#include <iostream>
#include <fstream>
#include <cstdlib>

#include "image.h"
#include "GradientLaplacian.h"

int main(int argc, char *argv[]) {
	int N = 256, M = 256, Q = 255;
	
	ImageType lenna(N, M, Q);
	readImage("../images/lenna.pgm", lenna);
	getPrewitt("lenna", lenna);
	//getSobel("lenna", lenna);
	//getLaplacian("lenna", lenna);
	
	ImageType sf(N, M, Q);
	readImage("../images/sf.pgm", sf);
	getPrewitt("sf", sf);
	//getSobel("lenna", sf);
	//getLaplacian("lenna", sf);

	return 0;
}
