#include <iostream>
#include <fstream>
#include <cstdlib>

#include "image.h"
#include "ImageSampling.h"

int main(int argc, char *argv[]) {

	std::cout << "Main called..." << std::endl;

	int N = 256;
	int M = 256;
	int Q = 255;
	int L = 256;

	// instantiate image objects
	ImageType lenna(N, M, Q);
	ImageType peppers(N, M, Q);

	// subsample original image
	imageSampling("lenna", lenna, 2);								// 256: factor of 2 to 128
	imageSampling("lenna", lenna, 4);								// 256: factor of 4 to 64
	imageSampling("lenna", lenna, 8);								// 256: factor of 8 to 32
	imageSampling("peppers", peppers, 2);
	imageSampling("peppers", peppers, 4);
	imageSampling("peppers", peppers, 8);
	// resize subsampled images
	imageSizing("lenna_ss2", lenna, 2);							// 128: factor of 2 back to 256
	imageSizing("lenna_ss4", lenna, 4);							// 64: factor of 4 back to 256
	imageSizing("lenna_ss8", lenna, 8);							// 32: factor of 8 back to 256
	imageSizing("peppers_ss2", peppers, 2);
	imageSizing("peppers_ss4", peppers, 4);
	imageSizing("peppers_ss8", peppers, 8);

	return 0;
}
