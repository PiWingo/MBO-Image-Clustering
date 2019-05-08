
#include "MBO_Object.h" 
#include "ggwFunc.h"
#include <cmath>

using namespace cv;
//using namespace std;

Mat histograma(Mat& imagem) {
	Mat hist;
	const float range1[] = { 0,256 };
	const float* range[] = { range1 };
	int histSize = 256;
	int channels[] = { 0 };
	calcHist(&imagem, 1, channels, Mat(), hist, 1, &histSize, range);
	for (int s = 0; s < 256; s++) {
		hist.at<float>(0, s) /= imagem.rows * imagem.cols;
	}
	return hist;
}

float entropiaMax(int estados) {
	float s = 0;
	for (int i = 0; i < estados; i++) {
		s += (1. / estados) * log(1. / estados);
	}
	return -s;
}

float entropia(Mat hist, int lim1, int lim2) {
	//normaliza
	float s = 0;
	for (int i = lim1; i < lim2; i++) {
		s += hist.at<float>(0, i);
	}
	if (s > 0)
		for (int i = lim1; i < lim2; i++) {
			hist.at<float>(0, i) /= s;
		}

	s = 0;
	for (int i = lim1; i < lim2; i++) {
		if (hist.at<float>(0, i) > 0)
			s += hist.at<float>(0, i) * log(hist.at<float>(0, i));
	}
	return -s;
}

float avaliacao(Mat & histograma, vector<int> cortes) {
	sort(cortes.begin(), cortes.end());

	float s = entropia(histograma, 0, cortes[0]);

	for (int i = 1; i < cortes.size(); i++) {
		s += entropia(histograma, cortes[i - 1], cortes[i]);
	}

	s += entropia(histograma, cortes.back(), 256);

	return s;
}

MBO_Object::MBO_Object(void)
{

}

MBO_Object::MBO_Object(int popsize, int dimension, int maxt)
	:Popsize(popsize), dim(dimension), max_t(maxt)
{

	

	period = 1.2; //com 1 e 0.5 embaixo n convergiu
	partition = 5.0 / 12; //isso aqui muda muito como converge
	BAR = partition;
	keep = 2;
	numButterfly1 = ceil(partition*Popsize);  // NP1 in paper
	numButterfly2 = Popsize - numButterfly1; // NP2 in paper
	maxStepSize = 1.0;
	MaxFEs = 1E4;

	ag.reserve(Popsize);
	ag1.reserve(numButterfly1);
	ag2.reserve(numButterfly2);
	Land1.reserve(dim);
	Land2.reserve(dim);
	Land1.assign(dim, 0.0);
	Land2.assign(dim, 0.0);
	Agent agent;
	for (int i = 0; i < Popsize; ++i)
	{
		agent.initialize(dim);
		ag.push_back(agent);
		if (i < numButterfly1)
		{
			ag1.push_back(agent);
		}
		else
		{
			ag2.push_back(agent);
		}
	}

}

MBO_Object::~MBO_Object(void)
{

}

void MBO_Object::MBO()
{
	vector<Agent> tempElitism;
	Agent bestAgent;

	// Write final results to output.txt
	ofstream fout;  //输出文件
	fout.open("output.txt", ios::out | ios::app);  //输出文件

	if (!fout.is_open())
	{
		cout << "Error opening file"; exit(1);
	}

	// 重定向
	streambuf *coutbackup;
	coutbackup = cout.rdbuf(fout.rdbuf());  //用 rdbuf() 重新定向

	Mat img = imread("Resources\\lenna.png", IMREAD_GRAYSCALE);
	Mat hist = histograma(img);

	CostFunction(ag, hist);
	PopSort(ag);
	bestAgent = ag[0];
	total_fit = ag[0].fit;


	// Begin the optimization loop
	for (int t = 0; t < max_t; ++t)
	{
		// Elitism strategy
		for (int ik = 0; ik < keep; ++ik)
		{
			tempElitism.push_back(ag.at(ik));
		}

		

		//////////////////    Divide the whole population into two subpopulations % % % %%%
		//	 Divide the whole population into Population1(Land1) and Population2(Land2)
		//	according to their fitness.
		//	 The monarch butterflies in Population1 are better than or equal to Population2.
		//	 Of course, we can randomly divide the whole population into Population1 and Population2.
		//	 We do not test the different performance between two ways.
		for (int popindex = 0; popindex < Popsize; ++popindex)
		{
			if (popindex < numButterfly1)
			{
				ag1[popindex] = ag[popindex];
			}
			else
			{
				ag2[popindex - numButterfly1] = ag[popindex];
			}
		}

		// Migration operator
		for (int k1 = 0; k1 < numButterfly1; ++k1)
		{
			for (int parnum1 = 0; parnum1 < dim; ++parnum1)
			{
				if (random(0.0, 1.0) * period <= partition)
				{
					int r2 = round(random(0, numButterfly1 - 1));
					Land1[parnum1] = ag1[r2].pos[parnum1];
				}
				else
				{
					int r3 = round(random(0, numButterfly2 - 1));
					Land1[parnum1] = ag2[r3].pos[parnum1];
				}
			}
			ag[k1].pos = Land1;
		}

		// Butterfly adjusting operator
		for (int k2 = 0; k2 < numButterfly2; ++k2)
		{
			auto scale = maxStepSize / pow(t + 1, 2); // Smaller step for local walk
			// R = exprnd(MU) returns an array of random numbers chosen from the
			// exponential distribution with mean parameter MU
			int StepSzie = ceil(exprnd(2 * max_t));
			auto delataX = LevyFlight(StepSzie, dim);

			for (int parnum2 = 0; parnum2 < dim; ++parnum2)
			{
				if (random(0.0, 1.0) >= partition)
				{
					Land2[parnum2] = bestAgent.pos[parnum2];
				}
				else
				{
					int r4 = round(random(0, numButterfly2 - 1));
					Land2[parnum2] = ag2[r4].pos[1];
					if (random(0.0, 1.0) > BAR)
					{
						Land2[parnum2] = Land2[parnum2] + scale * (delataX[parnum2] - 0.5);
					}
				}
			}
			ag[numButterfly1 + k2].pos = Land2;
		}

		FeasibleFunction(ag);
		CostFunction(ag,hist);
		PopSort(ag);
		bestAgent = ag[0];

		if (ag[0].fit < total_fit)
		{
			total_fit = ag[0].fit;
		}
		std::cout << total_fit << "   ";

		// Replace the worst with the previous generation's elites.
		for (int ik = 0; ik < keep; ++ik)
		{
			ag.at(Popsize - ik - 1) = tempElitism.at(ik);
		}
		tempElitism.clear();

	}  // end for max_t


	// Close the file
	fout.close();

	cout.rdbuf(coutbackup);  //取消，恢复屏幕输出
}

vector<double>  MBO_Object::LevyFlight(int StepSzie, int dim)
{
	// Allocate matrix for solutions
	vector<double>  delataX(dim);
	vector<double>  fx(StepSzie);

	// Loop over each dimension
	for (int i = 0; i < dim; ++i)
	{
		generate_n(fx.begin(), fx.size(),
			[=]()
		{
			return tan(pi*random(0.0, 1.0));
		});
		delataX[i] = sum(fx);
	}

	return delataX;
}

//calculate the fit according to pos-vector<double>
double MBO_Object::get_fit(vector<double> &pos)
{
	// Ackley Function
	double value, sum1 = 0.0, sum2 = 0.0;
	unsigned int i;
	double value1, value2, value3;

	for (i = 0; i < pos.size(); i++)
		sum1 += pow(pos[i], 2);

	for (i = 0; i < pos.size(); i++)
		sum2 += cos(2 * pi * pos[i]);

	value1 = -20 * exp(-0.2* sqrt(sum1 / pos.size()));
	value2 = -exp(sum2 / pos.size());
	value3 = 20 + exp(1); 

	value = value1 + value2 + value3;

	return value;
}

double MBO_Object::getEntropy(vector<double>& pos, Mat& image)
{
	int corteValor = (int)pos.at(0);
	vector<int> corte;
	corte.push_back(corteValor);
	double value = (double)avaliacao(image, corte);

	return value;
}


void MBO_Object::CostFunction(vector<Agent>& ag, Mat& image)
{
	for (int i = 0; i < Popsize; ++i)
	{
		ag.at(i).fit = getEntropy(ag.at(i).pos, image);
	}
}

//Compute the cost of each member in Population
void MBO_Object::FeasibleFunction(vector<Agent>& ag)
{
	for (int i = 0; i < Popsize; ++i)
	{
		for (int j = 0; j < dim; ++j)
		{
			ag.at(i).pos[j] = max(ag[i].pos[j], down);
			ag.at(i).pos[j] = min(ag[i].pos[j], up);
		}
	}
}

//Sort the population members from best to worst
void MBO_Object::PopSort(vector<Agent>& ag)
{	
	vector<double> cost(Popsize, 0.0);
	for (int i = 0; i < Popsize; ++i)
	{
		cost.at(i) = ag.at(i).fit;
	} 
	static int generation = 0;
	multimap<double, int> mm;
	for (unsigned int i = 0; i < cost.size(); ++i)
		mm.insert(make_pair(cost[i], i));

	vector<Agent>  agent(ag);

	auto cp = mm.begin();
	for (int i = 0; i < Popsize; ++i)
	{
		ag[i] = agent[(cp++)->second];
	}

	ofstream fout;  //output file
	fout.open("generation.txt", std::ios_base::app);  //??????

	if (!fout.is_open())
	{
		cout << "Error opening file"; exit(1);
	}
	streambuf* coutbackup;
	coutbackup = cout.rdbuf(fout.rdbuf());  //?? rdbuf() ???¶???

	for (int i = 0; i < ag.size(); i++) {
		cout << "ag fit: " << ag[i].fit;
		cout << " ag pos: " << ag[i].pos[1] <<endl;
	}
	cout << "-------fim da geracao: " << generation << endl;
	fout.close();

	cout.rdbuf(coutbackup);  //??????????????
	generation++;
}

double MBO_Object::sum(vector<double> & vec1)
{
	double ggwSum = 0.0;
	for (int ix = 0; ix != vec1.size(); ix++)
		ggwSum = ggwSum + vec1.at(ix);
	return ggwSum;
}

void MBO_Object::output()
{

	// Write final results to output.txt
	ofstream fout;  //输出文件
	fout.open("output.txt", ios::out);  //输出文件
	//fout.open("output.txt", ios::out | ios::app);  //输出文件

	if (!fout.is_open())
	{
		cout << "Error opening file"; exit(1);
	}

	// 重定向
	streambuf *coutbackup;
	coutbackup = cout.rdbuf(fout.rdbuf());  //用 rdbuf() 重新定向

	cout << endl << endl << "The final results are: " << endl;
	cout << "the parameter is : " << endl;
	cout << "popsize=" << Popsize << ", ";
	cout << "dimension=" << dim << ", ";
	cout << "maximum generation= " << max_t << endl;
	cout << "the result is: total_fit ---" << total_fit << endl;
	cout << "the solution is :" << endl;
	//for (int i = 0; i < dim; ++i)
	//{
	//	cout << ag[1].pos[i] << "   ";
	//}
	for_each(ag[1].pos.begin(), ag[1].pos.end(), [&](double a) { cout << a << " "; });

	cout << endl;

	// Close the file
	fout.close();

	cout.rdbuf(coutbackup);  //取消，恢复屏幕输出
}

