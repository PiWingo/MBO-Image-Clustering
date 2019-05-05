// main.cpp : Defines the entry point for the console application.
//

#include "MBO_Object.h"
#include <cstdlib>  
#include <iostream>
#include <ctime>
#include <iomanip>
#include <opencv2/opencv.hpp>

using namespace std;
using namespace cv;

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

void imhist(string const& name, Mat1b const& image) {
	// Set histogram bins count
	int bins = 256;
	int histSize[] = { bins };
	// Set ranges for histogram bins
	float lranges[] = { 0, 256 };
	const float* ranges[] = { lranges };
	// create matrix for histogram
	Mat hist;
	int channels[] = { 0 };

	// create matrix for histogram visualization
	int const hist_height = 256;
	Mat3b hist_image = Mat3b::zeros(hist_height, bins);

	calcHist(&image, 1, channels, Mat(), hist, 1, histSize, ranges, true, false);

	double max_val = 0;
	minMaxLoc(hist, 0, &max_val);

	// visualize each bin
	for (int b = 0; b < bins; b++) {
		float const binVal = hist.at<float>(b);
		int   const height = cvRound(binVal * hist_height / max_val);
		line
		(hist_image
			, Point(b, hist_height - height), Point(b, hist_height)
			, Scalar::all(255)
		);
	}
	imshow(name, hist_image);
}

int main()
{

	Mat image;
	image = imread("Resources\\pedagio.png", IMREAD_GRAYSCALE);

	if (image.empty()) {

	}
	else {
		Mat img;
		resize(image, img, Size(1024, 768), 0, 0, INTER_CUBIC);
		imshow("img hist", img);// Show our image inside it.
		imhist("image1 hist", img);
	}

	float melhor = 0;
	int itera = 0;

	
	clock_t start2 = clock();
	for (int i = 0; i < 256; i++) {

		vector<int> cortesHistograma = { i }; //a solução é descrita por esse vetor;
		Mat hist = histograma(image);

		float r = avaliacao(hist, cortesHistograma);
		//cout << "avaliacao: " << r << endl;
		if (melhor < r) {
			melhor = r;
			itera = i;
		}
	}
	clock_t end2 = clock();
	cout << "melhor avaliacao: " << melhor << endl;
	cout << "valor da faixa: " << itera << endl;
	cout << "the time is:  " << (double)(end2 - start2) / CLOCKS_PER_SEC << endl;

	srand((unsigned)time(NULL));

	int popsize = 0, maxt = 0, dim;
	// cin>>popsize>>dim>>n>>maxt>>cr>>f;
	popsize = 50;
	dim = 3;
	maxt = 50;

	MBO_Object mbo(popsize, dim, maxt);
	clock_t start = clock();
	mbo.MBO();
	clock_t end = clock();
	mbo.output();

	ofstream fout;  //output file
	fout.open("output.txt", ios::out | ios::app);  //??????

	if (!fout.is_open())
	{
		cout << "Error opening file"; exit(1);
	}

	streambuf* coutbackup;
	coutbackup = cout.rdbuf(fout.rdbuf());  //?? rdbuf() ???¶???

	cout << "CLOCKS_PER_SEC  " << CLOCKS_PER_SEC << endl;
	cout << "the time is:  " << (double)(end - start) / CLOCKS_PER_SEC << endl;

	// Close the file
	fout.close();

	cout.rdbuf(coutbackup);  //??????????????

	waitKey(0);
	return 0;
}

