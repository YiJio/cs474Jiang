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
	
	// apply blur and noise
	blurImage("lenna", 1);
	blurImage("lenna", 10);
	blurImage("lenna", 100);
	
	// apply inverse filtering with radius 25, 50, 100
	unblurImage("lenna_degraded1", 0, 25, 0);
	unblurImage("lenna_degraded1", 0, 50, 0);
	unblurImage("lenna_degraded1", 0, 100, 0);
	unblurImage("lenna_degraded10", 0, 25, 0);
	unblurImage("lenna_degraded10", 0, 50, 0);
	unblurImage("lenna_degraded10", 0, 100, 0);
	unblurImage("lenna_degraded100", 0, 25, 0);
	unblurImage("lenna_degraded100", 0, 50, 0);
	unblurImage("lenna_degraded100", 0, 100, 0);
	
	// apply wiener filtering with k 0.1, 0.01, 0.001 / 0.25, 0.025, 0.0025
	unblurImage("lenna_degraded1", 1, 0, 0.1);
	unblurImage("lenna_degraded1", 1, 0, 0.01);
	unblurImage("lenna_degraded1", 1, 0, 0.001);
	unblurImage("lenna_degraded10", 1, 0, 0.1);
	unblurImage("lenna_degraded10", 1, 0, 0.01);
	unblurImage("lenna_degraded10", 1, 0, 0.001);
	unblurImage("lenna_degraded100", 1, 0, 0.25);
	unblurImage("lenna_degraded100", 1, 0, 0.025);
	unblurImage("lenna_degraded100", 1, 0, 0.0025);
	
	return 0;
}
