#include<opencv2/opencv.hpp>
#include<iostream>

using namespace std;
using namespace cv;

void imHistogram(string const& name, Mat1b const& image) {
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
		int   const height = cvRound(binVal*hist_height / max_val);
		line
		(hist_image
			, Point(b, hist_height - height), Point(b, hist_height)
			, Scalar::all(255)
		);
	}
	imshow(name, hist_image);
}

/*int main() {
	string imageName;
	cout << "This simple code will convert any colored image to grayscale and output its histogram" << endl;
	cout << "Feel free to put any image file in the 'x64\\Release\\Resources\\' folder" << endl;
	cout << "You may use the 'image.jpg' file for testing" << endl;
	cout << "Type your file name (Example: image.jpg): ";
	cin >> imageName;
	cout << endl;
	Mat image1 = imread("Resources\\"+imageName, IMREAD_COLOR); // This will read a coloured image

	while (image1.empty()) {
		cout << "Error: Invalid File." << endl;
		cout << "Type your file name (Example: image.jpg): ";
		cin >> imageName;
		cout << endl;
		image1 = imread("Resources\\"+imageName, IMREAD_COLOR);
	}
	Mat1b image1_gray;
	Mat img;
	cvtColor(image1, image1_gray, COLOR_BGR2GRAY); // This will convert image1 to grayscale in image1_gray
	resize(image1_gray, img, Size(1024, 768), 0, 0, INTER_CUBIC); // Resizing the image1_gray to fit the screen
	imshow("image1",img );
	imhist("image1 hist",img);

	waitKey(0);


	return 0;

}*/