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
	srand((unsigned)time(NULL));

	int popsize = 0, maxt = 0, dim;
	// cin>>popsize>>dim>>n>>maxt>>cr>>f;
	popsize = 50;
	dim = 1;
	maxt = 10;

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

	return 0;
}

