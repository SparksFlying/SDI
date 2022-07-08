#pragma once
#include<iostream>
#include<iomanip>
#include<algorithm>
#include<cmath>
#include<random>
#include<numeric>
#include"segment.hpp"
#include"interval-tree/interval_tree.hpp"
using namespace std;
using namespace lib_interval_tree;
#define BLACK 0
#define RED 1
#define FORWARD 0
#define BACKWARD 1
using KT=uint64_t;
using PT=size_t;
const KT minval=(std::numeric_limits<KT>::min)();
const KT maxval=(std::numeric_limits<KT>::max)();
struct my_interval_t {
	KT low;
	KT high;
	segment<KT,PT,double>* p2seg;
};


int getHeight(interval_tree_t<KT>::node_type* root) {
		if (!root)
			return 0;
		int left_h = getHeight(root->left_);
		int right_h = getHeight(root->right_);
		return (left_h > right_h) ? (1 + left_h) : (1 + right_h);
}

// 将K个段填充为包含整个区间[minv,maxv]的intervals
// 根据flag为前向还是后向，在interval中保存右/左手边最靠近的段指针
vector<my_interval_t> padding_segments(std::vector<segment<KT,PT>>& segs,bool flag=FORWARD) {
		vector<my_interval_t> intervals;
		// 补充0~segs[0].start-1
		{
			if(flag==FORWARD){
				if(segs[0].start==0){

				}else{
					intervals.push_back(my_interval_t{0,segs[0].start-1,&segs[0]});
				}
			}else{
				if(segs[0].start==0){

				}else{
					intervals.push_back(my_interval_t{0,segs[0].start-1,new segment<KT,PT>(0,0,0,0,0,0)});
				}
			}
		}
		for (auto i = 0; i < segs.size(); ++i) {
			intervals.push_back(my_interval_t{ segs[i].start, segs[i].end ,&segs[i]});
			if (i < segs.size() - 1) {
				if (segs[i].end + 1 <= segs[i + 1].start - 1) {
					if(flag==FORWARD)
						intervals.push_back(my_interval_t{ segs[i].end + 1, segs[i + 1].start - 1,&segs[i+1] });
					else
						intervals.push_back(my_interval_t{ segs[i].end + 1, segs[i + 1].start - 1,&segs[i] });
				}
			}
		}
		// 补充segs.back().end+1~maxv
		{
			if(flag==FORWARD){
				intervals.push_back(my_interval_t{segs.back().end+1,std::numeric_limits<uint64_t>::max(),new segment<KT,PT>(0,0,0,0,std::numeric_limits<uint64_t>::max(),0)});
			}else{
				intervals.push_back(my_interval_t{segs.back().end+1,std::numeric_limits<uint64_t>::max(),&segs.back()});
			}
		}
		return intervals;
}



// 将线段树转为满二叉树，并返回数组
vector<interval_tree_t<KT>::node_type*> padding_intervaltree(interval_tree_t<KT>& tree) {
	using node_t=interval_tree_t<KT>::node_type;
	int height = getHeight(tree.root_);
	vector<interval_tree_t<KT>::node_type*> nodevec(std::pow(2, height) - 1, nullptr);
	nodevec[0] = tree.root_;
	for (auto i = 0; 2 * i + 2 < nodevec.size(); ++i) {
		if (nodevec[i]!=nullptr) {
			nodevec[2 * i + 1] = nodevec[i]->left_;
			nodevec[2 * i + 2] = nodevec[i]->right_;
		}
	}

	// 扩充为满二叉树的同时保留结点间关系
	KT minv = 1;
	constexpr KT maxv = (std::numeric_limits<KT>::max)();
	std::default_random_engine e;
	std::uniform_int_distribution<KT> u(minv, maxv);
	std::uniform_real_distribution<double> ud(0,1);
	for (auto i = 1; i < nodevec.size(); ++i) {
		if (nodevec[i] == nullptr) {
			// 随即生成伪结点
			size_t rand_number = u(e);
			nodevec[i] = new node_t(nodevec[size_t((i - 1) / 2)],node_t::interval_type{rand_number,rand_number - 1});
			nodevec[i]->data_ptr = new segment<>(rand_number,rand_number-1,0,0,0,0);
			nodevec[i]->left_=nullptr;nodevec[i]->right_=nullptr;
		}
	}
	for (auto i = 0; i < nodevec.size(); ++i) {
		if(2*i+1<nodevec.size()) nodevec[i]->left_=nodevec[2*i+1];
		if(2*i+2<nodevec.size()) nodevec[i]->right_=nodevec[2*i+2];
	}
	return nodevec;
}

// 仅平铺树结点，并返回数组
vector<interval_tree_t<KT>::node_type*> flatten_intervaltree(interval_tree_t<KT>& tree) {
	using node_t=interval_tree_t<KT>::node_type;
	int height = getHeight(tree.root_);
	vector<interval_tree_t<KT>::node_type*> nodevec(std::pow(2, height) - 1, nullptr);
	nodevec[0] = tree.root_;
	for (auto i = 0; 2 * i + 2 < nodevec.size(); ++i) {
		if (nodevec[i]!=nullptr) {
			nodevec[2 * i + 1] = nodevec[i]->left_;
			nodevec[2 * i + 2] = nodevec[i]->right_;
		}
	}
	return nodevec;
}