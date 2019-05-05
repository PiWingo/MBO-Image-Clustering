
#include "ggwFunc.h"
#include <iomanip>
#include <random>

// 生成1个随机数，范围(start, finish)，类型为double
double random(double start, double finish)
{
	return start + (finish - start)*rand() / (RAND_MAX + 1.0);
}

// 生成n个随机数，范围(start, finish)，类型为vector<double> 
vector<double>  random(double start, double finish, int n)
{
	vector<double>  vec(n);
	vec.reserve(n);
	auto ss = vec.size();
	for (unsigned int i = 0; i < vec.size(); ++i)
	{
		vec.at(i) = random(start, finish);
	}
	return vec;
}

// 生成n个随机数，范围(start, finish)，类型为double
/*
double exprnd(double lambda)
{
	std::random_device rd;
	std::mt19937 gen(rd());

	// if particles decay once per second on average,
	// how much time, in seconds, until the next one?
	std::exponential_distribution<> d(1.0/lambda);

	return d(gen);

}
*/

double exprnd(double lambda)
{
	double pV = 0.0;
	while (true)
	{
		pV = (double)rand() / (double)RAND_MAX;
		if (pV != 1)
		{
			break;
		}
	}
	//pV = (-1.0 / lambda)*log(1 - pV);
	pV = (-lambda)*log(1 - pV);
	return pV;
}
