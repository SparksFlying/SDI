#include"interface.hpp"
#include"protocols.hpp"
#include"config.hpp"
#include<vector>
#include<ctime>
#include<queue>
#include<functional>
using std::vector;
using std::queue;
using ophelib::Ciphertext;
using ophelib::PaillierFast;

#define LEFT 0
#define RIGHT 1

struct s_segment_node{
    Ciphertext low;
    Ciphertext high;
    Ciphertext max;
    s_segment_node* left;
    s_segment_node* right;
    rb_color color;
    std::function<Ciphertext(const ophelib::PaillierFast& crypto,const Ciphertext& input)> predictor;
    s_segment_node(){}
};

struct s_segment_tree{
    s_segment_node* root;
    int height;
};

// 加密线段树(线段树已被填充为满二叉树),并返回加密后的树
s_segment_tree encrypt_intervaltree(const ophelib::PaillierFast& crypto,vector<interval_tree_t<KT>::node_type*>& nodes){
    if(info){
            for(int i=0;i<nodes.size();++i){
                printf("\nori %d:[%u,%u,%u]\n",i,nodes[i]->low(),nodes[i]->high(),nodes[i]->max());
            }
    }
    
    vector<s_segment_node*> ctree(nodes.size());

    for(auto i=0;i<ctree.size();++i){
        ctree[i]=new s_segment_node;
        ctree[i]->left=nullptr;
        ctree[i]->right=nullptr;
    }

    // 加密数据以及保留结点间关系
    //#pragma omp parallel for
    for(auto i=0;i<nodes.size();++i){
        ctree[i]->low=crypto.encrypt(Integer(uint64_t(nodes[i]->low())));
        ctree[i]->high=crypto.encrypt(Integer(uint64_t(nodes[i]->high())));
        ctree[i]->max=crypto.encrypt(Integer(uint64_t(nodes[i]->max())));

        if(info){
            printf("\n%s->%u\n",Integer(uint64_t(nodes[i]->max())).to_string(false).c_str(),
                Integer(uint64_t(nodes[i]->max())).to_ulong());
            printf("\nenc %d:[%s,%s,%s]\n",i,Integer(uint64_t(nodes[i]->low())).to_string(false).c_str(),
                Integer(uint64_t(nodes[i]->high())).to_string(false).c_str(),
                Integer(uint64_t(nodes[i]->max())).to_string(false).c_str());    
        }
        if(info){
            printf("\ndec %d:[%s,%s,%s]\n",i,crypto.decrypt(ctree[i]->low).to_string(false).c_str(),
                crypto.decrypt(ctree[i]->high).to_string(false).c_str(),
                crypto.decrypt(ctree[i]->max).to_string(false).c_str());    
        }

        ctree[i]->color=nodes[i]->color();
        //ctree[i]->predictor
        segment<>* sptr=static_cast<segment<>*>(nodes[i]->data_ptr);
        if(sptr){
            Integer sl(uint64_t(pow(10,config::FLOAT_EXP)));
            sl*=(sptr->low_slope+sptr->high_slope)/2;
            Integer b(uint64_t(pow(10,config::FLOAT_EXP)));
            b*=sptr->pos;
            Integer tmp(uint64_t(pow(10,config::FLOAT_EXP)));
            tmp*=(sptr->low_slope+sptr->high_slope)/2;
            tmp*=sptr->start;
            b-=tmp;

            Ciphertext e_sl=crypto.encrypt(sl);
            Ciphertext e_b=crypto.encrypt(b);
            ctree[i]->predictor=[e_sl,e_b](const ophelib::PaillierFast& crypto,const Ciphertext& input)->Ciphertext{
                Ciphertext output=SM(crypto,e_sl,input).data*e_b.data;
                if(track) printf("\n%s * %s + %s = %s\n",crypto.decrypt(e_sl).to_string(false).c_str(),
                crypto.decrypt(input).to_string(false).c_str(),
                crypto.decrypt(e_b).to_string(false).c_str(),
                crypto.decrypt(output).to_string(false).c_str());
                return output;
            };
            if(info){
                printf("\n[sl=%s,b=%s]",sl.to_string(false).c_str(),
                        b.to_string(false).c_str());
                ctree[i]->predictor(crypto,crypto.encrypt(1));
            }
            
        }else{
            Ciphertext e_sl=crypto.encrypt(0);
            Ciphertext e_b=crypto.encrypt(0);
            ctree[i]->predictor=[e_sl,e_b](const ophelib::PaillierFast& crypto,const Ciphertext& input)->Ciphertext{
                Ciphertext output=SM(crypto,e_sl,input).data*e_b.data;
                return output;
            };
            if(info){
                printf("\n[sl=%s,b=%s]",Integer(0).to_string(false).c_str(),
                        Integer(0).to_string(false).c_str());
                ctree[i]->predictor(crypto,crypto.encrypt(1));
            }
        }
        if(2*i+1<nodes.size()) ctree[i]->left=ctree[2*i+1];
        if(2*i+2<nodes.size()) ctree[i]->right=ctree[2*i+2];
    }
    
    s_segment_tree res{ctree[0],static_cast<int>(std::log2(nodes.size()+1))};
    
    return res;
}

// 加密线段树(线段树为非满二叉树),并返回加密后的树
s_segment_tree encrypt_intervaltree_without_padding(const ophelib::PaillierFast& crypto,vector<interval_tree_t<KT>::node_type*>& nodes){
    Integer N=crypto.get_pub().n;
    vector<s_segment_node*> ctree(nodes.size());

    for(auto i=0;i<nodes.size();++i){
        if(nodes[i]!=nullptr) ctree[i]=new s_segment_node();
        else ctree[i]=nullptr;
    }
    // 加密数据以及保留结点间关系
    #pragma omp parallel for
    for(auto i=0;i<nodes.size();++i){
        if(nodes[i]!=nullptr)
        {
            ctree[i]->low=crypto.encrypt(Integer(nodes[i]->low()));
            ctree[i]->high=crypto.encrypt(Integer(nodes[i]->high()));
            ctree[i]->max=crypto.encrypt(Integer(nodes[i]->max()));
            ctree[i]->color=nodes[i]->color();
            //ctree[i]->predictor
            segment<>* sptr=static_cast<segment<>*>(nodes[i]->data_ptr);
            if(sptr){
                Integer sl(uint64_t(pow(10,config::FLOAT_EXP)));
                sl*=(sptr->low_slope+sptr->high_slope)/2;
                Integer b(uint64_t(pow(10,config::FLOAT_EXP)));
                b*=sptr->pos;
                Integer tmp(uint64_t(pow(10,config::FLOAT_EXP)));
                tmp*=(sptr->low_slope+sptr->high_slope)/2;
                tmp*=sptr->start;
                b-=tmp;

                Ciphertext e_sl=crypto.encrypt(sl);
                Ciphertext e_b=crypto.encrypt(b);
                ctree[i]->predictor=[e_sl,e_b](const ophelib::PaillierFast& crypto,const Ciphertext& input)->Ciphertext{
                    Ciphertext output=SM(crypto,e_sl,input).data*e_b.data;
                    if(track) printf("\n%s * %s + %s = %s\n",crypto.decrypt(e_sl).to_string(false).c_str(),
                    crypto.decrypt(input).to_string(false).c_str(),
                    crypto.decrypt(e_b).to_string(false).c_str(),
                    crypto.decrypt(output).to_string(false).c_str());
                    return output;
            };
            if(info){
                printf("\n[sl=%s,b=%s]",sl.to_string(false).c_str(),
                        b.to_string(false).c_str());
                ctree[i]->predictor(crypto,crypto.encrypt(1));
            }
            
            }else{
                Ciphertext e_sl=crypto.encrypt(0);
                Ciphertext e_b=crypto.encrypt(0);
                ctree[i]->predictor=[e_sl,e_b](const ophelib::PaillierFast& crypto,const Ciphertext& input)->Ciphertext{
                    Ciphertext output=SM(crypto,e_sl,input).data*e_b.data;
                    return output;
                };
                if(info){
                    printf("\n[sl=%s,b=%s]",Integer(0).to_string(false).c_str(),
                            Integer(0).to_string(false).c_str());
                    ctree[i]->predictor(crypto,crypto.encrypt(1));
                }
            }
            if(2*i+1<ctree.size()) ctree[i]->left=ctree[2*i+1];
            if(2*i+2<ctree.size()) ctree[i]->right=ctree[2*i+2];
        }
        
    }
    
    s_segment_tree res{ctree[0],getHeight(nodes[0])};
    return res;
}



// void inorder(PaillierFast& crypto,sIntervalTNode* root,vector<interval>& v){
//     if(!root) return;
//     Integer tmp=crypto.decrypt(root->inte.low.data%*crypto.get_n2());
//     Integer tmp2=crypto.decrypt(root->inte.high.data%*crypto.get_n2());
//     v.push_back(interval{crypto.decrypt(root->inte.low.data%(*crypto.get_n2())).to_ulong(),crypto.decrypt(root->inte.high.data%(*crypto.get_n2())).to_ulong(),nullptr});
//     inorder(crypto,root->left,v);
//     inorder(crypto,root->right,v);
// }
// void inorder(PaillierFast& crypto,IntervalTNode* root,vector<interval>& v){
//     if(!root) return;
//     v.push_back(interval{root->inte.low,root->inte.high,nullptr});
//     inorder(crypto,root->left,v);
//     inorder(crypto,root->right,v);
// }
// void verify(PaillierFast& crypto,sIntervalTree& ctree,IntervalTNode* root){
    
//     vector<interval> v1;
//     vector<interval> v2;
//     inorder(crypto,ctree.root,v1);
//     inorder(crypto,root,v2);
//     for(auto i=0;i<v1.size();++i){
//         if(v1[i].low==v2[i].low&&v1[i].high==v2[i].high){

//         }else{
//             printf("error!\n");
//         }
//     }
// }

Ciphertext SOverlap(const PaillierFast& crypto,const Ciphertext& a_low,const Ciphertext& a_high,const Ciphertext& b_low,const Ciphertext& b_high){
    Integer N=crypto.get_pub().n;
    // tmp1=E(a.high-b.low)
    Ciphertext tmp1=SMinus(crypto,a_high,b_low);
    Ciphertext part1=SLT(crypto,tmp1);
    // tmp2=E(b.high-a.low)
    Ciphertext tmp2=SMinus(crypto,b_high,a_low);
    Ciphertext part2=SLT(crypto,tmp2);
    return SMinus(crypto,crypto.encrypt(1),SOR(crypto,part1,part2));
    
}

bool SPathFlagTest(const PaillierFast& crypto,const Ciphertext& F,s_segment_node* node,const Ciphertext& input){
    Integer N=crypto.get_pub().n;
    if(!node->left){
        return RIGHT;
    }

    Ciphertext pf=SLT(crypto,node->left->max,input);
    pf=SMinus(crypto,crypto.encrypt(1),pf);
    
    int r=rand()%2;
    Ciphertext f=SM(crypto,SMinus(crypto,crypto.encrypt(1),F),pf).data*SM(crypto,F,crypto.encrypt(Integer(r))).data;

    return crypto.decrypt(f)==0?RIGHT:LEFT;
}

Ciphertext SOverlapTest(const PaillierFast& crypto,s_segment_node* node,const Ciphertext& input){
    if(info) printf("\nnode inte:%s,%s\n",
            crypto.decrypt(node->low).to_string(false).c_str(),
            crypto.decrypt(node->high).to_string(false).c_str());
    Ciphertext res=SOverlap(crypto,node->low,node->high,input,input);
    return res;
}

// 
Ciphertext Search(const PaillierFast& crypto,const s_segment_tree& tree,const Ciphertext& input){
    // // 记录从根到叶子沿途访问的结点
    // NTL::Vec<Ciphertext> path;
    // path.SetLength(tree.height);

    // // 记录沿途访问结点与i是否重叠
    // NTL::Vec<Ciphertext> overlapflags;
    // overlapflags.SetLength(tree.height);
    // long seq=0;
    // vector<int> plain_flags;
    // vector<int> debug_flags;
    Ciphertext P=crypto.encrypt(0);
    Ciphertext F=crypto.encrypt(0);

    // 遍历直到叶子结点
    s_segment_node* cur=tree.root;
    Integer N2=*crypto.get_n2();
    
    //printf("\n input:%s\n",crypto.decrypt(input.data).to_string(false).c_str());
    while(cur){
        Ciphertext f=SOverlapTest(crypto,cur,input);
        if(track){
            printf("\ncur=[%lu,%lu],input=%lu\n",crypto.decrypt(cur->low).to_ulong(),
                crypto.decrypt(cur->high).to_ulong(),
                crypto.decrypt(input).to_ulong());
        }
        
        F.data=F.data*f.data;
        
        Ciphertext acc=cur->predictor(crypto,input);
        P.data=P.data*SM(crypto,f,acc).data;

        if(track){
            printf("\nf=%d,F=%d,P=%lu\n",crypto.decrypt(f).to_int(),
                crypto.decrypt(F).to_int(),
                crypto.decrypt(P).to_ulong());
        }
        if(SPathFlagTest(crypto,F,cur,input)==LEFT){
            cur=cur->left;
        }else{
            cur=cur->right;
        }
        /*{
            printf("addr(cur)=%lu\n",crypto.decrypt(path[seq]).to_ulong());
            printf("i=[%lu,%lu],node range=[%lu,%lu]\n",crypto.decrypt(i.low).to_ulong(),crypto.decrypt(i.high).to_ulong(),crypto.decrypt(cur->inte.low).to_ulong(),crypto.decrypt(cur->inte.high).to_ulong());
            plain_flags.push_back(crypto.decrypt(overlapflags[seq]).to_int());
            if(crypto.decrypt(i.low).to_ulong()>=crypto.decrypt(cur->inte.low).to_ulong()&&
               crypto.decrypt(i.low).to_ulong()<=crypto.decrypt(cur->inte.high).to_ulong()
            ) debug_flags.push_back(1);
            else debug_flags.push_back(0);
        }*/
    }

    /* {
        printf("search done!\n");
        printf("palin flags:[");
        for(auto f:plain_flags){
            printf("%d,",f);
        }
        printf("]\n");
    }
    {
        printf("debug flags:[");
        for(auto f:debug_flags){
            printf("%d,",f);
        }
        printf("]\n");
    } */

    
    /* {
        Integer N2=*crypto.get_n2();
        auto start=crypto.decrypt(res->inte.p2seg->params[0].data%N2);
        auto end=crypto.decrypt(res->inte.p2seg->params[1].data%N2);
        auto sl=crypto.decrypt(res->inte.p2seg->params[2].data%N2);
        auto sp=crypto.decrypt(res->inte.p2seg->params[3].data%N2);
        auto pos=crypto.decrypt(res->inte.p2seg->params[4].data%N2);
        printf("\n");
    } */
    return P;
}


void delete_intervaltree(interval_tree_t<KT>::node_type* root){
    if(root){
        delete_intervaltree(root->left_);
        delete_intervaltree(root->right_);
        delete root;
    }
}

void delete_s_segment_tree(s_segment_node* root){
    if(root){
        delete_s_segment_tree(root->left);
        delete_s_segment_tree(root->right);
        delete root;
    }
}