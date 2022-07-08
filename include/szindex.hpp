#pragma once
#include<iostream>
#include"zindex.hpp"
#include"protocols.hpp"
#include"interface.hpp"
#include"ssegmentree.hpp"
#include<future>
#include<functional>
#define PARAL
#define DEBUG
#define PACK 0
#define NUM_THREADS 8

#ifdef DEBUG
#define DBGprint(...) printf(__VA_ARGS__)
#else
#define DBGprint(...)
#endif

#ifdef DEBUG
 TimerClock TC;
#endif

#ifdef DEBUG
#define TCstart \
    TC.update()
#else
#define TCstart
#endif

#ifdef DEBUG
#define TCend(...) \
    DBGprint(__VA_ARGS__,TC.getTimerMilliSec())
#else
#define TCend(...)
#endif

//基于已有的zindex加密后的安全索引
template<class KT,class PT,int dim,int bitperdim>
class SZIndex{
public:
    
    SZIndex(const ophelib::PaillierFast& crypto,ZIndex<KT,PT,dim,bitperdim>& index){
        this->N=crypto.get_pub().n;
        this->N2=*crypto.get_n2();
        this->epsilon=index.epsilon;

        // 记录zaddr的大小范围
        this->minz=crypto.encrypt(Integer(uint64_t(index.zdata[0])));
        this->maxz=crypto.encrypt(Integer(uint64_t(index.zdata.back())));

        // 记录原始数据集D
        //printf("encrypt data...\n");
        this->edata.SetDims(index.data.size(),dim);

        #pragma omp parallel for
        for(auto j=0;j<dim;++j){
            for(auto i=0;i<index.data.size();++i){
                this->edata[i][j]=crypto.encrypt(index.data[i][j]);
            }
        }
        //printf("encrypt data done!\n");

        //printf("build tree...\n");
        TCstart;
        
        {
            #pragma omp parallel for
            for(int i=FORWARD;i<=BACKWARD;++i)
            {
                if(i==FORWARD) build_tree(crypto,index,this->FTree,FORWARD,true);
                else build_tree(crypto,index,this->BTree,BACKWARD,true);
            }
            
        }
        // 构造前向加密线段树
        //build_tree(crypto,index,this->FTree,FORWARD);
        // 构造后向加密线段树
        //build_tree(crypto,index,this->BTree,BACKWARD);
        TCend("build tree time:%f ms\n");
    }
    void build_tree(const PaillierFast& crypto,ZIndex<KT,PT,dim,bitperdim>& index,s_segment_tree& ctree,bool flag=FORWARD,bool padding_flag=true){
        lib_interval_tree::interval_tree_t<uint64_t> tree;
//#ifdef DEBUG_
        //printf("\tpadding segments...\n");
        //printf("\t\tinit size of segs is %lu\n",index.segs.size());
//#endif
        vector<my_interval_t> interval_vec=padding_segments(index.segs,flag);
//#ifdef DEBUG_
        //printf("\t\tafter padding, size of intervals is %lu\n",intervals.size());
//#endif
        /* {
            size_t count=0;
            for(auto& seg:index.segs){
                printf("seg %lu=[%lu,%lu]\n",count++,seg.start,seg.end);
            }
            count=0;
            for(auto& inte:intervals){
                printf("interval %lu=[%lu,%lu]\n",count++,inte.low,inte.high);
            }
        } */
        //printf("\tinserting intervals...\n");
        for(auto& inte: interval_vec){
            tree.insert({inte.low,inte.high},(void*)inte.p2seg);
        }
        //printf("\theight of tree is %d\n",getHeight(tree.root));
        /* {
            auto res=IntervalT_Search(&tree,interval{2550816767065905,2550816767065905,nullptr});
            printf("res find [%u,%u]\n",res->inte.low,res->inte.high);
        } */
        //printf("\tpadding interval tree...\n");
        vector<interval_tree_t<uint64_t>::node_type*> all_nodes;
        if(padding_flag) all_nodes=padding_intervaltree(tree);
        else all_nodes=flatten_intervaltree(tree);

        //printf("\tencrypt interval tree...\n");
        if(info){
            for(int i=0;i<all_nodes.size();++i){
                printf("\nori %d:[%u,%u,%u]\n",i,all_nodes[i]->low(),all_nodes[i]->high(),all_nodes[i]->max());
            }
        }
        if(padding_flag) ctree=encrypt_intervaltree(crypto,all_nodes);
        else ctree=encrypt_intervaltree_without_padding(crypto,all_nodes);

        /* {
            verify(crypto,ctree,all_nodes[0]);
            printf("encrypt done!\n");
        } */
            
    }

    // std::pair<PT,PT> point_query(ophelib::PaillierFast& crypto,ophelib::Ciphertext ekey){
    //     // 计时器
    //      TC;
    //     // 存储ekey与各段的所属关系
    //     NTL::Vec<ophelib::Ciphertext> cmp_res;
    //     cmp_res.SetLength(esegs.NumRows());
    //     Integer N=crypto.get_pub().n;
    //     Integer N2=*crypto.get_n2();
    //     // 检测ekey与各段的所属关系，属于则为1，否则为0
    //     TCstart
    //     {
    //         Random& r=Random::instance();
    //         //#pragma omp parallel for
    //         for(auto i=0;i<esegs.NumRows();++i){
    //             Integer ra=r.rand_int(2);//掷一个随机数属于{0,1}
    //             // part1=E(2Q+1-2start)
    //             Integer part1=ekey.data*
    //                           ekey.data*
    //                           crypto.encrypt(1).data*
    //                           SM(crypto,crypto.encrypt(ophelib::Integer(-1)%N),esegs[i][start].data*esegs[i][start].data).data;
    //             // part2=E(2Q-2end-1)
    //             Integer part2=ekey.data*
    //                           ekey.data*
    //                           SM(crypto,crypto.encrypt(ophelib::Integer(-1)%N),esegs[i][end].data*esegs[i][end].data).data*
    //                           crypto.encrypt(Integer(-1)%N).data;
    //             // tmp=E[(2Q+1-2start)(2Q-2end-1)]
    //             Integer tmp=SM(crypto,part1,part2).data;

    //             // 随机翻转
    //             if(ra==1){
    //                 tmp=SM(crypto,crypto.encrypt(Integer(-1)%N),tmp).data;
    //             }
    //             cmp_res[i]=SLT(crypto,tmp);
    //             if(ra==1){
    //                 cmp_res[i]=crypto.encrypt(1).data*SM(crypto,crypto.encrypt(Integer(-1)%N),cmp_res[i]).data;
    //             }
    //         }
    //     }
    //     auto detect_segs_time=TC.getTimerMilliSec();
    //     /* std::cout<<"[";
    //     for(auto i=0;i<esegs.NumRows();++i){
    //         std::cout<<crypto.decrypt(cmp_res[i].data%N2)<<",";
    //     }
    //     std::cout<<"]"<<std::endl; */
        
    //     // 获取ekey所对应的加密段
    //     TCstart
    //     NTL::Vec<ophelib::Ciphertext> eseg=dot(crypto,ophelib::Vector::transpose(esegs),cmp_res);
    //     auto get_eseg_time=TC.getTimerMilliSec();
    //     //在段内寻找

    //     TCstart
    //     //ls_base=E(pos*10^lsp)
    //     //hs_base=E(pos*10^hsp)
    //     ophelib::Ciphertext ls_base=SM(crypto,eseg[pos],crypto.encrypt(Integer(10).pow(crypto.decrypt(eseg[lsp].data%N2).to_long())%N));
    //     ophelib::Ciphertext hs_base=SM(crypto,eseg[pos],crypto.encrypt(Integer(10).pow(crypto.decrypt(eseg[hsp].data%N2).to_long())%N));
    //     //predict_pos=slope*(key-start)+pos
    //     ophelib::Integer pred_pos;
    //     {
    //         ophelib::Integer lp=crypto.decrypt
    //                             (
    //                                 ls_base.data*
    //                                 SM(
    //                                     crypto,
    //                                     eseg[ls],
    //                                     ekey.data*SM(crypto,crypto.encrypt(Integer(-1)%N),eseg[start]).data
    //                                 ).data%N2
    //                             )/Integer(10).pow(crypto.decrypt(eseg[lsp].data%N2).to_long());
    //         ophelib::Integer hp=crypto.decrypt
    //                             (
    //                                 hs_base.data*
    //                                 SM(
    //                                     crypto,
    //                                     eseg[hs],
    //                                     ekey.data*SM(crypto,crypto.encrypt(Integer(-1)%N),eseg[start]).data
    //                                 ).data%N2
    //                             )/Integer(10).pow(crypto.decrypt(eseg[hsp].data%N2).to_long());
    //         pred_pos=(lp+hp)/2;
    //     }
    //     Integer lp=pred_pos-crypto.decrypt(epsilon.data%N2);
    //     Integer hp=pred_pos+crypto.decrypt(epsilon.data%N2);
    //     auto pred_pos_time=TC.getTimerMilliSec();
    //     lp=lp<0?0:lp;
    //     hp=hp.to_long()>this->edata.NumRows()?Integer(this->edata.NumRows()-1):hp;
    //     //std::cout<<std::endl;
    //     /* std::cout<<"low pos:"<<lp<<std::endl
    //              <<"high pos:"<<hp<<std::endl; */

    //     // 抽取[lp,hp]之间的edata
    //     std::cout<<"detect segs relationship time: "<<detect_segs_time<<" ms"<<std::endl;
    //     std::cout<<"get eseg time: "<<get_eseg_time<<" ms"<<std::endl;
    //     std::cout<<"predict pos time: "<<pred_pos_time<<" ms"<<std::endl;
        
    //     return std::make_pair<PT,PT>(lp.to_int(),
    //                                  hp.to_int());
    // }
    std::vector<PT> range_query(const ophelib::PaillierFast& crypto,const std::vector<std::vector<Ciphertext>>& eRange,const EQueryRectangle<Ciphertext>& eQR,uint64_t param1=0,uint64_t param2=0){
        vector<PT> res;
        for(auto subrange:eRange){
            //printf("\n input:%s,%s\n",crypto.decrypt(subrange[0]),crypto.decrypt(subrange[1]));
            //对查询范围作缩放,两个端点必定属于一个interval
            /* {
                printf("before shringking QR=[%lu,%lu]\n",crypto.decrypt(subQR[0].data%N2).to_ulong(),crypto.decrypt(subQR[1].data%N2).to_ulong());
            } */
            // subQR[0]=SMAX(crypto,subQR[0],this->minz);
            // subQR[1]=SMIN(crypto,subQR[1],this->maxz);
            /* {
                printf("after shringking QR=[%lu,%lu]\n",crypto.decrypt(subQR[0].data%N2).to_ulong(),crypto.decrypt(subQR[1].data%N2).to_ulong());
                printf("shringk range done!\n");
            } */
            Ciphertext e_start,e_end;
            TCstart;
#ifdef PARAL
            {
                //printf("parallal!\n");
                auto func1=std::bind(Search,crypto,this->FTree,subrange[0]);
                auto f1=std::async(std::launch::async, func1);
                auto func2=std::bind(Search,crypto,this->BTree,subrange[1]);
                auto f2=std::async(std::launch::async, func2);
                e_start=f1.get();
                e_end=f2.get();
            }
#else
            {
                e_start=Search(crypto,this->FTree,subrange[0]);
                e_end=Search(crypto,this->BTree,subrange[1]);
            }
#endif

            {
                PT start;
                {
                    Integer tmp=crypto.decrypt(e_start);
                    tmp/=uint64_t(pow(10,config::FLOAT_EXP));
                    tmp-=this->epsilon;
                    if(tmp<0) tmp=0;
                    start=tmp.to_ulong();
                }
                
                PT end;
                {
                    Integer tmp=crypto.decrypt(e_end);
                    tmp/=uint64_t(pow(10,config::FLOAT_EXP));
                    tmp+=this->epsilon;
                    end=tmp.to_ulong();
                    end=min(end,PT(this->edata.NumRows()-1));
                }
                
                {
                    start=param1;
                    end=param2;
                }
                if(0){
                    printf("\nszindex predicts:[%u,%u]\n",start,end);
                }

                TCend("first stage cost=%f ms\n");
                if(dim==1){
                    for(PT i=start;i<=end;++i){
                        res.push_back(i);
                    }
                    return res;
                }
                
                if(PACK!=1){
                    TCstart;
#ifdef PARAL
                    #pragma omp parallel for //num_threads(NUM_THREADS)
#endif
                    for(PT i=start;i<=end;++i){
                        //TCstart
                        if(eQR.isFallin(crypto,this->edata[i])){
                            res.push_back(i);
                        }
                        //cout<<TC.getTimerMilliSec()<<endl;
                    }
                    TCend("second stage cost=%f ms\n");
                }
                else{
                    TCstart;
                    res=packing_search(start,end,crypto,eQR);
                    TCend("second stage cost=%f ms\n");
                }
            }
            //res.push_back((crypto.decrypt(start_pos.data%N2)/crypto.decrypt(start_seg->params[sp].data%N2)).to_uint());
            //res.push_back((crypto.decrypt(end_pos.data%N2)/crypto.decrypt(end_seg->params[sp].data%N2)).to_uint());
            return res;
        }
        
    }

    std::vector<PT> packing_search(const PT& start,const PT& end,const PaillierFast& crypto,const EQueryRectangle<Ciphertext>& eQR){
        if(end<start) return vector<PT>{};
        std::vector<Vec<Ciphertext>> P(dim);
        Vec<Vec<PackedCiphertext>> Packed_data;
        Vec<Vec<PackedCiphertext>> Pack_QR_low;
        Vec<Vec<PackedCiphertext>> Pack_QR_high;
        Vec<Vec<PackedCiphertext>> Packed_part;
        Pack_QR_low.SetLength(dim);
        Pack_QR_high.SetLength(dim);
        Packed_data.SetLength(dim);
        Packed_part.SetLength(dim);
        for(auto& p:P) p.SetLength(end-start+1);
        for(int i=0;i<dim;++i){
            Pack_QR_low[i]=Vector::pack_ciphertexts_vec(Vec<Ciphertext>(NTL::INIT_SIZE_TYPE{},end-start+1,eQR.get_minvec()[i]),64,crypto);
            Pack_QR_high[i]=Vector::pack_ciphertexts_vec(Vec<Ciphertext>(NTL::INIT_SIZE_TYPE{},end-start+1,eQR.get_maxvec()[i]),64,crypto);
        }

        for(int i=0;i<dim;++i){
            P[i].SetLength(end-start+1);
            for(size_t j=start;j<=end;++j){
                P[i][j-start]=edata[j][i];
            }
            //auto tmppp=Vector::pack_ciphertexts_vec(P[i],64, crypto);
            Packed_data[i]=Vector::pack_ciphertexts_vec(P[i],64, crypto);
            //Packed_data[i].move(tmppp);
        }

        //
        Vec<Integer> tmp_res(NTL::INIT_SIZE_TYPE{},end-start+1,0);
        //vector<Vec<Integer>> random_noise(dim,Vector::rand_bits_neg(end-start+1,8);
        #pragma omp parallel for
        for(int i=0;i<dim;++i){
            Packed_part[i].SetLength(Packed_data[i].length());
            size_t idx=0;
            for(size_t j=0;j<Packed_data[i].length();++j){
                PackedCiphertext tmp1(SMinus(crypto,Packed_data[i][j].data.data,Pack_QR_low[i][j].data.data),Packed_data[i][j].n_plaintexts,64);
                PackedCiphertext tmp2(SMinus(crypto,Pack_QR_high[i][j].data.data,Packed_data[i][j].data.data),Packed_data[i][j].n_plaintexts,64);
                Vec<Integer> tmp1_unpack=Vector::decrypt_pack(tmp1,crypto);
                Vec<Integer> tmp2_unpack=Vector::decrypt_pack(tmp2,crypto);
                // Vec<Integer> debug1_unpack=Vector::decrypt_pack(Packed_data[i][j],crypto);
                // Vec<Integer> debug2_unpack=Vector::decrypt_pack(Pack_QR_low[i][j],crypto);
                // Vec<Integer> debug3_unpack=Vector::decrypt_pack(Pack_QR_high[i][j],crypto);
                //cout<<"data:"<<debug1_unpack;
                // cout<<"\n------------------\n";
                // cout<<"QR low:"<<debug2_unpack;
                // cout<<"\n------------------\n";
                // cout<<"QR high:"<<debug3_unpack;
                // cout<<"\n------------------\n";
                // cout<<"data-low:"<<tmp1_unpack;
                // cout<<"\n------------------\n";
                // cout<<"high-data"<<tmp2_unpack;
                //cout<<endl;
                for(size_t k=0;k<tmp1_unpack.length();++k){
                    if(tmp1_unpack[k]>=0&&tmp2_unpack[k]>=0) tmp_res[idx+k]+=1;
                }
                idx+=tmp1_unpack.length();
            }
        }
        std::vector<PT> res;
        for(size_t k=0;k<end-start+1;++k){
            if(tmp_res[k]==dim) res.push_back(k+start);
        }

        return res;
    }
public:
    // 存储加密后的原始数据
    NTL::Mat<Ciphertext> edata;
    // 存储管理段集的加密线段树
    // 前向树
    s_segment_tree FTree;
    // 后向树
    s_segment_tree BTree;
    // 存储预测误差参数
    size_t epsilon;
    // 存储zdddr的范围端点值
    Ciphertext minz;
    Ciphertext maxz;

    // 
    size_t start=0; //start key
    size_t end=1;   //end key
    size_t sl=2;    //slope
    size_t sp=3;    //slope precision
    size_t pos=4;   //start pos
    size_t len=5;

    Integer N;
    Integer N2;
    bool isParal=false;
public:
    ~SZIndex(){
        //delete_sintervaltree(this->FTree.root);
        //delete_sintervaltree(this->BTree.root);
    }
};



