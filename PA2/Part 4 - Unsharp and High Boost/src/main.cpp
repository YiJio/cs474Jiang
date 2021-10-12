#include <iostream>
#include <fstream>
#include <cstdlib>

#include "image.h"
#include "UnsharpHighBoost.h"

int main(int argc, char *argv[]) {
	int N = 256, M = 256, Q = 255;
	
	// test lenna image
	ImageType lenna(N, M, Q);
	readImage("../images/lenna.pgm", lenna);
	getGaussian("lenna", lenna, 7);								// perform gaussian for smooth image
	ImageType lennaLP(N, M, Q);
	readImage("../images/lenna_gaussian_mask7.pgm", lennaLP);	// read and store smooth image
	getSharpMask("lenna", lenna, lennaLP);						// generate sharp mask
	ImageType lennaMask(N, M, Q);
	readImage("../images/lenna_sharp_mask.pgm", lennaMask);		// read and store sharp mask
	getSharp("lenna", lenna, lennaMask, 1);						// perform sharp filtering k=1
	getSharp("lenna", lenna, lennaMask, 2);						// perform sharp filtering k>1
	
	// test f_16 image
	ImageType f16(N, M, Q);
	readImage("../images/f_16.pgm", f16);
	getGaussian("f_16", f16, 7);
	ImageType f16LP(N, M, Q);
	readImage("../images/f_16_gaussian_mask7.pgm", f16LP);	
	getSharpMask("f_16", f16, f16LP);
	ImageType f16Mask(N, M, Q);
	readImage("../images/f_16_sharp_mask.pgm", f16Mask);
	getSharp("f_16", f16, f16Mask, 1);
	getSharp("f_16", f16, f16Mask, 2);

	return 0;
}
