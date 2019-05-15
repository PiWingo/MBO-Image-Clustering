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
#include <opencv2/opencv.hpp>
using namespace std;
using namespace cv;

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
	double getEntropy(double pos, const Mat& hist);
	void CostFunction(vector<Agent>& ag, const Mat& hist);
	void PopSort(vector<Agent>& ag);
	void FeasibleFunction(vector<Agent>& ag);
	vector<double> LevyFlight(int StepSzie, int dim);
	double sum(vector<double> & vec1);
	void output();
	MBO_Object(void);
	MBO_Object(int num, int dimension, int maxt);
public:
	~MBO_Object(void);
};
