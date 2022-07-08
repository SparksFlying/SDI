#pragma once
#include<vector>
#include"segment.hpp"
#include"utility.hpp"
#include<string>
#include<libzinc/zinc.hh>
#include<cstdio>
#include<cstdint>
#include<cassert>
#include<cmath>
#include<fstream>
#include<numeric>
#include"segs_rmi.h"
#include<algorithm>
#include<functional>
#include<unordered_map>
using std::string;
using std::pair;
using std::vector;
using std::ofstream;
const size_t NOT_FOUND=0xFFFFFFFF;



// 2-6 dim Z-curve Index
template<class KT=uint32_t,class PT=size_t,int dim=2,int bitperdim=32>
class ZIndex{
public:
    ZIndex(std::vector<std::vector<KT>>& data,PT epsilon=32);
    //将段信息写入文件用于训练RMI
    void write(const string& path);
    //item属于segs[0].start到segs[M-1].end之间，寻找item所属的段标范围
    std::pair<PT,PT> find_segment_range(const vector<KT>& item);
    //寻找item所属的段,若item不属于任何一个段，则返回距离它最近的右边的一个段
    PT find_segment(const vector<KT>& item,bool& isFound);
    //寻找item在data中所属的范围
    std::pair<PT,PT> find_range(PT seg_pos,uint64_t key);
    //若item属于某个段，则寻找item在data中的位置(序号)；否则，返回位置的最大下界
    PT point_query(const vector<KT>& item);

    std::vector<std::vector<KT>> range_query(const QueryRectangle<KT>& QR);
    std::vector<PT> range_query_v2(const QueryRectangle<KT>& QR);
    //threshold表示启用跳跃机制的最大错误阈值
    std::vector<PT> range_query_v2_1(const QueryRectangle<KT>& QR,short threshold=10);
    std::vector<PT> range_query_v3(const QueryRectangle<KT>& QR);

    std::vector<std::vector<KT>> sequential_search(const PT& start,const PT& end);
    
public:
    std::vector<segment<uint64_t,PT,double>> segs;
    std::vector<std::vector<KT>> data;
    std::vector<uint64_t> zdata;
    // 记录zdata中元素pos与data中对应元素pos之间的关系
    std::unordered_map<PT,PT> loc_table;
    PT epsilon;
};

template<class KT,class PT,int dim,int bitperdim>
ZIndex<KT,PT,dim,bitperdim>::ZIndex(std::vector<std::vector<KT>>& data,PT epsilon):epsilon(epsilon),zdata(data.size()),data(data.size(),vector<KT>(dim)){
    vector<pair<uint64_t,size_t>> tmp_arr(data.size());
    for(size_t i=0;i<data.size();++i){
        tmp_arr[i].first=morton_code<dim,bitperdim>::encode(vec2arr<KT,dim>(data[i]));
        tmp_arr[i].second=i;
    }
    std::sort(tmp_arr.begin(),tmp_arr.end(),[](const pair<uint64_t,size_t>& a,const pair<uint64_t,size_t>& b){
        return a.first<b.first;
    });

    for(size_t i=0;i<data.size();++i){
        this->zdata[i]=tmp_arr[i].first;
        for(int j=0;j<dim;++j) this->data[i][j]=data[tmp_arr[i].second][j];
    }
    
    std::vector<PT> pos(data.size());
    for(auto i=0;i<data.size();++i) pos[i]=i;

    //拟合线性段
    this->segs=shringkingCone(zdata,pos,epsilon);

    while(this->segs.size()>350){
        epsilon+=16;
        this->segs=shringkingCone(zdata,pos,epsilon);
    }

    //加载segs_rmi参数
    segs_rmi::load("/home/zanglang/source/SecureLearnedIndex/RMI/rmi_data");

    // 
}

template<class KT,class PT,int dim,int bitperdim>
void ZIndex<KT,PT,dim,bitperdim>::write(const string& path){
    ofstream fs(path,std::ios::binary);
    
    if(!fs.is_open()){
        std::cout<<"can't open file!"<<std::endl;
        exit(-1);
    }
    uint64_t size=segs.size();
    fs.write((char*)&size,sizeof(size));
    for(auto& seg:segs){
        fs.write((char*)&seg.start,sizeof(seg.start));
    }
}

template<class KT,class PT,int dim,int bitperdim>
std::pair<PT,PT> ZIndex<KT,PT,dim,bitperdim>::find_segment_range(const vector<KT>& item){
    uint64_t key=morton_code<dim,bitperdim>::encode(vec2arr<KT,dim>(item));

    //通过RMI预测位置
    size_t err;
    size_t pos=segs_rmi::lookup(key,&err);
    size_t low=pos>err?pos-err:0;
    size_t high=pos+err>=segs.size()?segs.size()-1:pos+err;
    return std::pair<PT,PT>{low,high};
}

template<class KT,class PT,int dim,int bitperdim>
PT ZIndex<KT,PT,dim,bitperdim>::find_segment(const vector<KT>& item,bool& isFound){
    uint64_t key=morton_code<dim,bitperdim>::encode(vec2arr<KT,dim>(item));
    
    //case 1
    if(key<segs[0].start){
        isFound=false;
        return 0;
    }
    //case 2，由于右边没有段，因此返回溢出
    if(key>segs.back().end){
        isFound=false;
        return segs.size();
    }
    //说明segs[0].start<=key<=segs[M-1].end,那么一定可以找到一个key所属的最大下界段
    //那么剩下仅有两种情况，1是key落在某个段内，对应case 3，2是key落在段与段的间隙中，对应case 4
    //case 3
    //s[i].start<=key<=s[i].end then return s[i]
    //case 4
    //s[i].end<key<s[i+1].start then return s[i+1]
    std::pair<PT,PT> range=find_segment_range(item);
    for(PT i=range.first;i<=range.second;++i){
        if(key>=segs[i].start&&key<=segs[i].end){
            isFound=true;
            return i;
        }
        if(i+1<=range.second&&key>segs[i].end&&key<segs[i+1].start){
            isFound=false;
            return i+1;
        }
    }
    return NOT_FOUND;
}

template<class KT,class PT,int dim,int bitperdim>
std::pair<PT,PT> ZIndex<KT,PT,dim,bitperdim>::find_range(PT seg_pos,uint64_t key){
    PT pred_pos=segs[seg_pos].pos+
                    (key-segs[seg_pos].start)*
                    (segs[seg_pos].high_slope+segs[seg_pos].low_slope)/2;
    PT  low=pred_pos>segs[seg_pos].pos+epsilon?
            pred_pos-epsilon:
            segs[seg_pos].pos;
    PT high=pred_pos+epsilon;
            high=high>=segs[seg_pos].pos+segs[seg_pos].len?
            segs[seg_pos].pos+segs[seg_pos].len-1:
            high;
    return std::make_pair<PT,PT>(std::move(low),std::move(high));
}

template<class KT,class PT,int dim,int bitperdim>
PT ZIndex<KT,PT,dim,bitperdim>::point_query(const vector<KT>& item){
    uint64_t key=morton_code<dim,bitperdim>::encode(vec2arr<KT,dim>(item));

    bool isFound;
    //TimerClock TC;
    //TC.update();
    PT seg_pos=find_segment(item,isFound);
    //std::cout<<"find "<<key<<"'s segment cost:"<<TC.getTimerMicroSec()<<"μs"<<std::endl;
    //item不属于任何一个段，应返回其位置的最大下界
    if(!isFound){
        //item超出zdata的最大值
        if(seg_pos==segs.size()){
            return data.size();
        }
        else{
            return segs[seg_pos].pos;
        }
    }
    //item处于segs[seg_pos]所表示的范围内，但是并不一定属于data
    else{
        //TC.update();
        auto pred_pos_range=find_range(seg_pos,key);
        //std::cout<<"find "<<key<<"'s range cost:"<<TC.getTimerMicroSec()<<"μs"<<std::endl;
        //TC.update();
        auto it=std::distance(zdata.begin(),std::lower_bound(zdata.begin()+pred_pos_range.first,zdata.begin()+pred_pos_range.second,key));
        //std::cout<<"find "<<key<<"'s iterator cost:"<<TC.getTimerMicroSec()<<"μs"<<std::endl;
        return it;
    }
    
}

//拆分QR为多个连续区间，实际运行时间过长，废弃
/* template<class KT,class PT>
std::vector<std::vector<KT>> ZIndex<KT,PT,dim,bitperdim>::range_query(const std::vector<std::vector<KT>>& QR){
    uint64_t start_key=morton_code<dim,bitperdim>::encode(std::array<KT,2>{QR[0][0],QR[0][1]});
    uint64_t end_key=morton_code<dim,bitperdim>::encode(std::array<KT,2>{QR[1][0],QR[1][1]});
    zinc::morton::AABB<2, 32> aabb {start_key, end_key};
    zinc::morton::region<2, 32> region = aabb.to_intervals();

    vector<vector<KT>> res;
    #pragma omp parallel
    for(auto& interval:region.intervals){
        auto start_p=morton_code<dim,bitperdim>::decode(interval.start);
        auto end_p=morton_code<dim,bitperdim>::decode(interval.end);
        auto start_range=point_query(vector<KT>{start_p[0],start_p[1]});
        auto end_range=point_query(vector<KT>{end_p[0],end_p[1]});
        for(auto i=start_range.first;i<end_range.second;++i){
            auto code=morton_code<dim,bitperdim>::encode(std::array<KT,2>{this->data[i][0],this->data[i][1]}).data;
            if(code>=interval.start&&code<=interval.end){
                res.push_back(this->data[i]);
            }
        }
    }
    return res;
} */

//预测QR起始点与终点的位置，遍历并测试
template<class KT,class PT,int dim,int bitperdim>
std::vector<PT> ZIndex<KT,PT,dim,bitperdim>::range_query_v2(const QueryRectangle<KT>& QR){
    uint64_t start_key=morton_code<dim,bitperdim>::encode(vec2arr<KT,dim>(QR.get_minvec()));
    uint64_t end_key=morton_code<dim,bitperdim>::encode(vec2arr<KT,dim>(QR.get_maxvec()));

    vector<PT> res;
    //TimerClock TC;

    //TC.update();
    auto start_pos=point_query(QR.get_minvec());
    //std::cout<<"find "<<start_key<<"'s start_pos "<<start_pos<<" cost:"<<TC.getTimerMilliSec()<<"ms"<<std::endl;
    //TC.update();
    //auto end_pos=point_query(std::vector<KT>{QR[0][1],QR[1][1]});
    //std::cout<<"find "<<end_key<<"'s end_pos "<<end_pos<<" cost:"<<TC.getTimerMilliSec()<<"ms"<<std::endl;
    //起始键已超出范围，返回空
    if(start_pos==data.size()){
        return res;
    }
    //std::cout<<"total search num: "<<end_pos-start_pos<<std::endl;
    /* std::cout<<"v2 first z addr is "<<start_key
             <<", last z addr is "<<end_key
             <<std::endl;
    std::cout<<"pred range:";
    std::cout<<start_pos<<","
             <<end_pos<<std::endl; */
    //TC.update();
    for(auto i=start_pos;i<data.size()&&zdata[i]<=end_key;++i){
        if(QR.isFallin(this->data[i])){
            res.push_back(i);
        }
    }
    //std::cout<<"loop cost:"<<TC.getTimerMilliSec()<<"ms"<<std::endl;
    return res;
}

//预测QR起始点与终点的位置，遍历并测试(加入跳跃机制)
template<class KT,class PT,int dim,int bitperdim>
std::vector<PT> ZIndex<KT,PT,dim,bitperdim>::range_query_v2_1(const QueryRectangle<KT>& QR,short threshold){
    uint64_t start_key=morton_code<dim,bitperdim>::encode(vec2arr<KT,dim>(QR.get_minvec()));
    uint64_t end_key=morton_code<dim,bitperdim>::encode(vec2arr<KT,dim>(QR.get_maxvec()));
    zinc::morton::AABB<dim,bitperdim> region{start_key,end_key};

    vector<PT> res;
    //TimerClock TC;
    //TC.update();
    auto start_pos=point_query(QR.get_minvec());
    //std::cout<<"find start_pos "<<start_pos<<" cost:"<<TC.getTimerMilliSec()<<"ms"<<std::endl;
    //TC.update();
    //auto end_pos=point_query(std::vector<KT>{QR[0][1],QR[1][1]});
    //std::cout<<"find end_pos "<<end_pos<<" cost:"<<TC.getTimerMilliSec()<<"ms"<<std::endl;
    //起始键已超出范围，返回空
    if(start_pos==data.size()){
        return res;
    }
    //std::cout<<"total search num: "<<end_pos-start_pos<<std::endl;
    /* std::cout<<"v2 first z addr is "<<start_key
             <<", last z addr is "<<end_key
             <<std::endl;
    std::cout<<"pred range:";
    std::cout<<start_pos<<","
             <<end_pos<<std::endl; */
    //TC.update();
    short count=0;
    for(auto i=start_pos;i<data.size()&&zdata[i]<=end_key;++i){
        if(QR.isFallin(this->data[i])){
            res.push_back(i);
        }else{
            count++;
            //跳跃到下一个位置
            if(count==threshold){
                auto next_addr=region.morton_get_next_address();
                while(zdata[i]<=next_addr.first){
                    if(QR.isFallin(this->data[i])){
                        res.push_back(i);
                    }
                    ++i;
                }
                while(zdata[i]<next_addr.second) ++i;
                // 在for循环自增前测试
                if(QR.isFallin(this->data[i])){
                        res.push_back(i);
                }
                count=0;
            }
        }
    }
    //std::cout<<"loop cost:"<<TC.getTimerMilliSec()<<"ms"<<std::endl;
    return res;
}

//原始循环
template<class KT,class PT,int dim,int bitperdim>
std::vector<PT> ZIndex<KT,PT,dim,bitperdim>::range_query_v3(const QueryRectangle<KT>& QR){
    vector<PT> res;
    //vector<PT> seqs;
    //#pragma omp parallel for
    for(auto i=0;i<this->data.size();++i){
        if(QR.isFallin(this->data[i])){
            res.push_back(i);
            //seqs.push_back(i);
        }
    }
    /* std::cout<<"loop based range:";
    std::cout<<seqs[0]<<","<<seqs.back(); */
    return res;
}

