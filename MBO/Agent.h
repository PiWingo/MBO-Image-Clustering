#pragma once

#include "ggwFunc.h" 
#include <cstdlib>
#include <vector>
#include <ctime>
using namespace std;

const double up = 256.00;
const double down = 0.00;

class Agent
{
private:
	int dim;
public:
	double fit;
	vector<double> pos;
	void initialize(int d);
public:
	Agent(void);
public:
	~Agent(void);
};
