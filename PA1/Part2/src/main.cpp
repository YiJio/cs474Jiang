#include <iostream>
#include <fstream>
#include <cstdlib>

#include "image.h"
#include "ImageQuantization.h"

int main(int argc, char *argv[]) {

	std::cout << "Main called..." << std::endl;

	int N = 256;
	int M = 256;
	int Q = 255;
	int L = 256;

	// instantiate image objects
	ImageType lenna(N, M, Q);
	ImageType peppers(N, M, Q);

	// quantize original image
	imageQuantization("lenna", lenna, 7, 0);				// bpp=7 -> L=128
	imageQuantization("lenna", lenna, 5, 0);				// bpp=5 -> L=32
	imageQuantization("lenna", lenna, 3, 0);				// bpp=3 -> L=8
	imageQuantization("lenna", lenna, 1, 0);				// bpp=1 -> L=2
	imageQuantization("peppers", peppers, 7, 0);			// bpp=7 -> L=128
	imageQuantization("peppers", peppers, 5, 0);			// bpp=5 -> L=32
	imageQuantization("peppers", peppers, 3, 0);			// bpp=3 -> L=8
	imageQuantization("peppers", peppers, 1, 0);			// bpp=1 -> L=2
	// subsample quantized images
	imageQuantization("lenna_L128", lenna, 7, 1);
	imageQuantization("lenna_L32", lenna, 5, 1);
	imageQuantization("lenna_L8", lenna, 3, 1);
	imageQuantization("lenna_L2", lenna, 1, 1);
	imageQuantization("peppers_L128", peppers, 7, 1);
	imageQuantization("peppers_L32", peppers, 5, 1);
	imageQuantization("peppers_L8", peppers, 3, 1);
	imageQuantization("peppers_L2", peppers, 1, 1);

	return 0;
}
