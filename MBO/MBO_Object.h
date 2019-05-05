#pragma once 

#include "Agent.h"
#include "ggwFunc.h"
#include <cmath>
#include <algorithm>
#include <vector>
#include <iostream>
#include <fstream>
#include <array>
#include <map>
#include<numeric>
#include<cassert>
#include <memory>
#include <stdio.h>
#include <stdlib.h>
#include <time.h> 
using namespace std;

class MBO_Object
{
private:
	vector<Agent> ag;
	vector<Agent> ag1;
	vector<Agent> ag2;
	int Popsize;
	int dim;
	double period;
	double partition;
	double BAR;
	int keep;
	int numButterfly1;  // NP1 in paper
	int numButterfly2; // NP2 in paper
	vector<double> Land1;
	vector<double> Land2;
	double maxStepSize;
	int max_t;
	int MaxFEs;
public:
	double total_fit;
	void MBO();
	double get_fit(vector<double> &pos);
	void CostFunction(vector<Agent>& ag);
	void PopSort(vector<Agent>& ag);
	void FeasibleFunction(vector<Agent>& ag);
	vector<double> LevyFlight(int StepSzie, int dim);
	double sum(vector<double> & vec1);
	void input();
	void output();
	MBO_Object(void);
	MBO_Object(int num, int dimension, int maxt);
public:
	~MBO_Object(void);
};
