#include <iostream>
#include <fstream>
#include <cstdlib>

#include "image.h"
#include "ImageSampling.h"
#include "ImageQuantization.h"

int main(int argc, char *argv[]) {

	std::cout << "Main called..." << std::endl;

	int N = 256;
	int M = 256;
	int Q = 255;
	ImageType img(N,M,Q);

	/* part 1 */
	// subsample original image
	imageSampling("../images/lenna", img, 2);							// 256: factor of 2 to 128
	imageSampling("../images/lenna", img, 4);							// 256: factor of 4 to 64
	imageSampling("../images/lenna", img, 8);							// 256: factor of 8 to 32
	imageSampling("../images/peppers", img, 2);
	imageSampling("../images/peppers", img, 4);
	imageSampling("../images/peppers", img, 8);
	// resize subsampled images
	imageSizing("../images/lenna_ss2", img, 2);						// 128: factor of 2 back to 256
	imageSizing("../images/lenna_ss4", img, 4);						// 64: factor of 4 back to 256
	imageSizing("../images/lenna_ss8", img, 8);						// 32: factor of 8 back to 256
	imageSizing("../images/peppers_ss2", img, 2);
	imageSizing("../images/peppers_ss4", img, 4);
	imageSizing("../images/peppers_ss8", img, 8);

	/* part 2 */
	// quantize original image
	imageQuantization("../images/lenna", img, 7, 0);			// bpp=7 -> L=128
	imageQuantization("../images/lenna", img, 5, 0);			// bpp=5 -> L=32
	imageQuantization("../images/lenna", img, 3, 0);			// bpp=3 -> L=8
	imageQuantization("../images/lenna", img, 1, 0);			// bpp=1 -> L=2
	imageQuantization("../images/peppers", img, 7, 0);			// bpp=7 -> L=128
	imageQuantization("../images/peppers", img, 5, 0);			// bpp=5 -> L=32
	imageQuantization("../images/peppers", img, 3, 0);			// bpp=3 -> L=8
	imageQuantization("../images/peppers", img, 1, 0);			// bpp=1 -> L=2
	// subsample quantized images
	imageQuantization("../images/lenna_L128", img, 7, 1);
	imageQuantization("../images/lenna_L32", img, 5, 1);
	imageQuantization("../images/lenna_L8", img, 3, 1);
	imageQuantization("../images/lenna_L2", img, 1, 1);
	imageQuantization("../images/peppers_L128", img, 7, 1);
	imageQuantization("../images/peppers_L32", img, 5, 1);
	imageQuantization("../images/peppers_L8", img, 3, 1);
	imageQuantization("../images/peppers_L2", img, 1, 1);

	return 0;
}
