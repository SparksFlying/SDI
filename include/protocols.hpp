#pragma once
#include "ophelib/paillier_fast.h"
#include "ophelib/vector.h"
#include "ophelib/omp_wrap.h"
#include "ophelib/packing.h"
#include "ophelib/util.h"
#include "ophelib/ml.h"
#include "ophelib/random.h"
#include <vector>
#include <libzinc/encoding.hh>
#include "utility.hpp"
using namespace ophelib;
using std::vector;
// ae=E(a),be=E(b)
// res= a op b
// 安全乘法协议
// E(res)=E(a*b)
Ciphertext SM(const PaillierFast& crypto,const Ciphertext& ae,const Ciphertext& be);

// 安全或协议
// E(res)=E(a|b)
Ciphertext SOR(const PaillierFast& crypto,const Ciphertext& ae,const Ciphertext& be);

// 安全异或协议
// E(res)=E(a^b)
Ciphertext SXOR(const PaillierFast& crypto,const Ciphertext& ae,const Ciphertext& be);

// 安全min协议
// E(res)=E(min(a,b))
Ciphertext SMIN(const PaillierFast& crypto,const Ciphertext& ae,const Ciphertext& be);

// 安全位分解协议
// 0<=a<2^m,return <E(a0),E(a1)...E(am-1)>
Vec<Integer> SBD(const PaillierFast& crypto, Ciphertext& ae,size_t m);

// 安全减法协议
// return E(a-b)
Ciphertext SMinus(const PaillierFast& crypto,const Ciphertext& ae,const Ciphertext& be);

// 安全翻转协议
// if ae=E(1),return E(0)
// if ae=E(0),return E(1)
Ciphertext SNOT(const PaillierFast& crypto,const Ciphertext& ae);

// 安全比较协议
// if a<0,return 1
// else return 0
Ciphertext SLT(const PaillierFast& crypto,const Ciphertext& ae);

// 矩阵-向量点积
NTL::Vec<ophelib::Ciphertext> dot(const ophelib::PaillierFast& crypto,const Mat<ophelib::Ciphertext>& A,const NTL::Vec<ophelib::Ciphertext>& B);

// 安全查询协议
// 查询点p是否在查询矩形Q内部
// 返回明文
template<int dim,int bitperdim>
bool SQuery(const ophelib::PaillierFast& crypto,const vector<Ciphertext>& R,const NTL::Vec<Ciphertext>& p){
    auto R_min_vec=morton_code<dim,bitperdim>::decode(crypto.decrypt(R[0].data%(*crypto.get_n2())).to_ulong());
    auto R_max_vec=morton_code<dim,bitperdim>::decode(crypto.decrypt(R[1].data%(*crypto.get_n2())).to_ulong());
    for(size_t d=0;d<dim;++d){
        uint32_t val=crypto.decrypt(p[d].data%(*crypto.get_n2())).to_uint();

        if(val<R_min_vec[d]||val>R_max_vec[d]){
            return false;
        }
    }
    return true;
}
//==========================================================================

//ae=E(a),be=E(b),res=E(a*b)
Ciphertext SM(const PaillierFast& crypto,const Ciphertext& ae,const Ciphertext& be){
    //global
    Integer N=crypto.get_pub().n;
    Integer N2=*crypto.get_n2();
    Ciphertext res;
    //return crypto.encrypt((crypto.decrypt(ae.data%(*crypto.get_n2()))*crypto.decrypt(be.data%(*crypto.get_n2())))%N);
    //P1 does
    Random& r=Random::instance();
    Integer ra=r.rand_int(N);
    Integer rb=r.rand_int(N);
    //std::string ra_str=ra.to_string_(),rb_str=rb.to_string_();
    Ciphertext tmp_a=crypto.encrypt(ra).data*ae.data;
    Ciphertext tmp_b=crypto.encrypt(rb).data*be.data;
    //P2 does
    Integer ha=crypto.decrypt(tmp_a.data%N2);
    Integer hb=crypto.decrypt(tmp_b.data%N2);
    Integer h=(ha*hb)%N;
    Ciphertext tmp_h=crypto.encrypt(h);
    //P1 does
    auto s=tmp_h.data*ae.data.pow_mod_n(N-rb,N2);
    auto tmp_s=s*be.data.pow_mod_n(N-ra,N2);
    res=tmp_s*crypto.encrypt(ra*rb%N).data.pow_mod_n(N-1,N2);
    /* std::cout<<"a="<<crypto.decrypt(ae.data%N2)<<", b="<<crypto.decrypt(be.data%N2)<<std::endl
             <<"a*b="<<crypto.decrypt(res.data%N2)<<std::endl; */
    return res;
}



//xe=E(x),0<=x<2^m,return <E(x0),E(x1)...E(xm-1)>
// 辅助协议
Integer EncLSB(const PaillierFast& crypto,Integer& T,size_t i);
bool SVR(const PaillierFast& crypto, Ciphertext& xe,Vec<Integer> xe_bits, size_t m=64);

Vec<Integer> SBD(const PaillierFast& crypto, Ciphertext& xe,size_t m){
    //global
    Integer N=crypto.get_pub().n;
    Integer N2=*crypto.get_n2();

    Vec<Integer> res;
    res.SetLength(m);
    Integer l=Integer(2).inv_mod_n(N);
    step:
    Integer T=xe.data;
    for(size_t i=0;i<m;++i){
        res[i]=EncLSB(crypto,T,i);
        Integer Z=(T*res[i].pow_mod_n(N-1,N2)).pow_mod_n(1,N2);
        T=Z.pow_mod_n(l,N2);
    }
    bool gama=SVR(crypto,xe,res,m);
    if(gama){
        return res;
    }else{
        goto step;
    }
}
bool SVR(const PaillierFast& crypto, Ciphertext& xe,Vec<Integer> xe_bits, size_t m){
    //global
    Integer N=crypto.get_pub().n;
    Integer N2=N.pow(2);
    Random& r=Random::instance();
    //P1 does
    Integer U=1;
    for(size_t i=0;i<m;++i){
        U*=xe_bits[i].pow_mod_n(Integer(2).pow(i),N2);
    }
    U=U%N2;
    Integer V=(U*xe.data.pow_mod_n(N-1,N2)).pow_mod_n(1,N2);
    Integer W=V.pow_mod_n(r.rand_int(N),N2);
    //P2 does
    if(crypto.decrypt(W%(*crypto.get_n2()))==0){
        return 1;
    }else{
        return 0;
    }
}
Integer EncLSB(const PaillierFast& crypto,Integer& T,size_t i){
    //global
    Integer N=crypto.get_pub().n;
    Integer N2=N.pow(2);
    Random& r=Random::instance();

    //P1 does
    Integer r1=r.rand_int(N);
    Integer Y=(T*crypto.encrypt(r1).data).pow_mod_n(1,N2);
    //P2 does
    Integer y=crypto.decrypt(Y%(*crypto.get_n2()));
    Integer alpha;
    if(y%Integer(2)==0){
        alpha=crypto.encrypt(0).data;
    }else{
        alpha=crypto.encrypt(1).data;
    }
    //P1 does
    if(r1%Integer(2)==0){
        return alpha;
    }else{
        return (crypto.encrypt(1).data*alpha.pow_mod_n(N-1,N2)).pow_mod_n(1,N2);
    }
}


//if a<b return ae,else return be
Ciphertext SMIN(const PaillierFast& crypto,const Ciphertext& ae,const Ciphertext& be){
    //global
    Integer N=crypto.get_pub().n;
    Integer N2=*crypto.get_n2();

    //P1 does
    Ciphertext tmp=ae.data*be.data.pow_mod_n(N-1,N2);
    //P2 does
    Ciphertext cmp=SLT(crypto,tmp);
    //P1 dose
    Integer p1=SM(crypto,ae,cmp).data;
    Integer p2=SM(crypto,be,SNOT(crypto,cmp)).data;
    //std::cout<<p1.to_string_()<<endl<<p2.to_string_()<<endl;
    Integer p=(p1*p2).pow_mod_n(1,N2);
    return p;

}

// max(a,b)=a+b-min(a,b)
Ciphertext SMAX(const PaillierFast& crypto,const Ciphertext& ae,const Ciphertext& be){
    return SMinus(crypto,ae.data*be.data,SMIN(crypto,ae,be));
    // if(crypto.decrypt(ae.data%(*crypto.get_n2()))>crypto.decrypt(be.data%(*crypto.get_n2()))){
    //     return ae;
    // }else{
    //     return be;
    // }
}


// Secure Less Than Protocol
// if res<0 return 1 otherwise return 0
// 未实现
Ciphertext SLT(const PaillierFast& crypto,const Ciphertext& res){
    Integer tmp=crypto.decrypt(res.data%(*crypto.get_n2()));
    if(tmp<0){
        return crypto.encrypt(1);
    }else{
        return crypto.encrypt(0);
    }
    // //P1 does
    // int r=rand()%2;
    // if(r==0){
    //     //P2 does
    //     Integer tmp=crypto.decrypt(res.data%(*crypto.get_n2()));
    // }else{
        
    // }
}

// Secure Less Than Or Equal Protocol
// if res<=0 return 1 otherwise return 0
// 未实现
Ciphertext SLTOE(const PaillierFast& crypto,const Ciphertext& res){
    Integer tmp=crypto.decrypt(res.data%(*crypto.get_n2()));
    if(tmp<=0){
        return crypto.encrypt(1);
    }else{
        return crypto.encrypt(0);
    }
}

Ciphertext SLT(const PaillierFast& crypto,const Ciphertext& a,const Ciphertext& b){
    Ciphertext tmp=SMinus(crypto,a,b);
    Integer res=crypto.decrypt(tmp.data%(*crypto.get_n2()));
    if(res<0){
        return crypto.encrypt(1);
    }else{
        return crypto.encrypt(0);
    }
}


Ciphertext SNOT(const PaillierFast& crypto,const Ciphertext& ae){
    return SMinus(crypto,crypto.encrypt(Integer(1)),ae);
}


Ciphertext SEQ(const PaillierFast& crypto,const Ciphertext& ae,const Ciphertext& be){
    return SXOR(crypto,SLT(crypto,SMinus(crypto,ae,be)),SLT(crypto,SMinus(crypto,be,ae)));
    //return crypto.encrypt(crypto.decrypt(a.data%(*crypto.get_n2()))==crypto.decrypt(b.data%(*crypto.get_n2())));
}

// res=E(a-b)
Ciphertext SMinus(const PaillierFast& crypto,const Ciphertext& a,const Ciphertext& b){
    Integer N2=*crypto.get_n2();
    return a.data*b.data.pow_mod_n(crypto.get_pub().n-1,N2);
}

// E(res)=E(a|b)
// res=a+b-a*b
Ciphertext SOR(const PaillierFast& crypto,const Ciphertext& ae,const Ciphertext& be){
    return SMinus(crypto,ae.data*be.data,SM(crypto,ae,be));
    // if(crypto.decrypt(ae.data%(*crypto.get_n2()))==0&&crypto.decrypt(be.data%(*crypto.get_n2()))==0){
    //     return crypto.encrypt(0);
    // }else
    // {
    //     return crypto.encrypt(1);
    // }
}

// E(res)=E(a ^ b)
// res=(1-a)*b+(1-b)*a
Ciphertext SXOR(const PaillierFast& crypto,const Ciphertext& ae,const Ciphertext& be){
    return SM(crypto,SMinus(crypto,crypto.encrypt(1),ae),be).data*SM(crypto,SMinus(crypto,crypto.encrypt(1),be),ae).data;
}


// 辅助协议
Ciphertext aux_dot(const ophelib::PaillierFast& crypto,const Vec<Ciphertext> &A, const Vec<Ciphertext> &B) {
            const long n = A.length();
            if (n != B.length())
                dimension_mismatch();
            if (n == 0)
                error_exit("empty vector");

            Ciphertext ret = SM(crypto,A[0],B[0]);
            #pragma omp parallel for
            for(long i = 1; i < n; i++) {
                ret.data *= SM(crypto,A[i],B[i]).data;
            }

            return ret;
        }
// 矩阵-向量点积
NTL::Vec<ophelib::Ciphertext> dot(const ophelib::PaillierFast& crypto,const Mat<ophelib::Ciphertext>& A,const NTL::Vec<ophelib::Ciphertext>& B){
    const long n = A.NumRows(),d = A.NumCols();

    if(d != B.length())
        dimension_mismatch();
    if(n == 0 || d == 0)
        error_exit("empty matrix");

    Vec<Ciphertext> ret;
    ret.SetLength(n);

    omp_set_nested(0);
    #pragma omp parallel for
    for(long i = 0; i < n; i++) {
        ret[i] = aux_dot(crypto,A[i], B);
    }

    return ret;
}


// QR -> eRange
template<int dim=2,int bitperdim=32>
std::vector<ophelib::Ciphertext> encrypt_rectangle(const ophelib::PaillierFast& crypto,const QueryRectangle<uint32_t>& QR){
    uint64_t start=morton_code<dim,bitperdim>::encode(vec2arr<uint32_t,dim>(QR.get_minvec()));
    uint64_t end=morton_code<dim,bitperdim>::encode(vec2arr<uint32_t,dim>(QR.get_maxvec()));
    Ciphertext a=crypto.encrypt(Integer(start)%crypto.get_pub().n);
    Ciphertext b=crypto.encrypt(Integer(end)%crypto.get_pub().n);
    std::vector<ophelib::Ciphertext> res{a,b};
    return res;
}



// 查询矩形,矩形格式为[xmin,ymin,zmin,xmax,ymax,zmax]
template<class KT=Ciphertext>
class EQueryRectangle{
public:
	EQueryRectangle(const vector<KT>& vals):dim(vals.size()/2),minvec(vector<KT>(dim)),maxvec(vector<KT>(dim)){
		for(auto i=0;i<dim;++i){
			minvec[i]=vals[i];
			maxvec[i]=vals[i+dim];
		}
	}
    template<class plaintext_type>
    EQueryRectangle(const PaillierFast& crypto,const QueryRectangle<plaintext_type>& QR):dim(QR.dim),minvec(vector<KT>(dim)),maxvec(vector<KT>(dim)){
        for(auto i=0;i<QR.dim;++i){
            minvec[i]=crypto.encrypt(QR.get_minvec()[i]);
            maxvec[i]=crypto.encrypt(QR.get_maxvec()[i]);
        }
    }
	const vector<KT>& get_minvec()const{
		return minvec;
	}
	const vector<KT>& get_maxvec()const{
		return maxvec;
	}
	bool isFallin(const PaillierFast& crypto,const NTL::Vec<Ciphertext>& p)const{
        Ciphertext res=crypto.encrypt(0);
		for(auto i=0;i<dim;++i){
			// if(crypto.decrypt(p[i])<crypto.decrypt(minvec[i])||crypto.decrypt(p[i])>crypto.decrypt(maxvec[i])){
			// 	return false;
			// }
            Integer tmp_cmp=SM(crypto,SMinus(crypto,p[i],minvec[i]),SMinus(crypto,p[i],maxvec[i])).data;
            Integer f=SLT(crypto,tmp_cmp.pow_mod_n(crypto.get_pub().n-1,*crypto.get_n2())).data;
            //Integer f=SOR(crypto,SLT(crypto,tmp_cmp),SEQ(crypto,tmp_cmp,crypto.encrypt(0))).data;
            res=res.data*SMinus(crypto,crypto.encrypt(1),f).data;
		}
		return crypto.decrypt(res)==dim?true:false;
	}
private:
	int dim;
	vector<KT> minvec;
	vector<KT> maxvec;
};