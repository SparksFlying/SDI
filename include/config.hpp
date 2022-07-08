#pragma once
#include<iostream>
#include<map>
using std::map;
#define info 0
#define track 0
namespace config{
    const size_t EPSILON=32;
    const uint64_t SIZE=1000;
    const int FLOAT_EXP=15;
    map<int,int> dim_bitperdim_table{
        {1,64},
        {2,32},
        {3,21},
        {4,16},
        {5,12},
        {6,10}
    };
};
