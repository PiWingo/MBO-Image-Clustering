
#pragma once

#include <random>
#include<algorithm>
#include<vector>
#include<cassert>
#include<iostream>
#include <type_traits>
using namespace std;

// the definition for some constants
const double pi = 3.14159265;


// 最小的数
const double epsilonDouble = numeric_limits<double>::epsilon(); // double型
const int epsilonInt = numeric_limits<int>::epsilon(); // int型


// 生成1个随机数，范围(start, finish)，类型为double
double random(double start, double finish); 

// 生成n个随机数，范围(start, finish)，类型为double
vector<double>  random(double start, double finish, int n);

// 生成n个指数分布随机数，类型为double
double exprnd(double lambda);
