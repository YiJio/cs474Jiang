#include <iostream>
#include <fstream>
#include <cstdlib>

#include "image.h"
#include "fft.h"

int main(int argc, char *argv[]) {
	int N = 256, M = 256, Q = 255;
	
	// generate lenna image to test
	ImageType lenna(N, M, Q);
	readImage("../images/lenna.pgm", lenna);
	
	// apply blur
	blurImage("lenna", lenna, 1);
	blurImage("lenna", lenna, 10);
	blurImage("lenna", lenna, 100);
	
	// read blurred images
	ImageType lenna1(N, M, Q);
	ImageType lenna10(N, M, Q);
	ImageType lenna100(N, M, Q);
	readImage("../images/lenna_blur1.pgm", lenna1);
	readImage("../images/lenna_blur10.pgm", lenna10);
	readImage("../images/lenna_blur100.pgm", lenna100);
	
	// apply inverse filtering with radius 5, 10, 15
	unblurImage("lenna_blur1", lenna1, 0, 10, 0);
	unblurImage("lenna_blur1", lenna1, 0, 50, 0);
	unblurImage("lenna_blur1", lenna1, 0, 100, 0);
	unblurImage("lenna_blur10", lenna10, 0, 10, 0);
	unblurImage("lenna_blur10", lenna10, 0, 50, 0);
	unblurImage("lenna_blur10", lenna10, 0, 100, 0);
	unblurImage("lenna_blur100", lenna100, 0, 10, 0);
	unblurImage("lenna_blur100", lenna100, 0, 50, 0);
	unblurImage("lenna_blur100", lenna100, 0, 100, 0);
	
	// apply wiener filtering with k 0.1, 0.01, 0.001 / 0.25, 0.025, 0.0025
	unblurImage("lenna_blur1", lenna1, 1, 0, 0.1);
	unblurImage("lenna_blur1", lenna1, 1, 0, 0.01);
	unblurImage("lenna_blur1", lenna1, 1, 0, 0.001);
	unblurImage("lenna_blur10", lenna10, 1, 0, 0.1);
	unblurImage("lenna_blur10", lenna10, 1, 0, 0.01);
	unblurImage("lenna_blur10", lenna10, 1, 0, 0.001);
	unblurImage("lenna_blur100", lenna100, 1, 0, 0.25);
	unblurImage("lenna_blur100", lenna100, 1, 0, 0.025);
	unblurImage("lenna_blur100", lenna100, 1, 0, 0.0025);
	
	return 0;
}
