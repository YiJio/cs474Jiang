#include <iostream>
#include <fstream>
#include <cstdlib>

#include "image.h"
#include "Correlation.h"

int main(int argc, char *argv[]) {
	int N = 442;
	int M = 288;
	int Q = 255;
	int Nm = 83;
	int Mm = 55;
	
	ImageType image(N, M, Q);
	ImageType mask(N, M, Q);	
	readImage("../images/Image.pgm", image);
	readImage("../images/Pattern.pgm", mask);
	getCorrelation("Image", image, mask);

	return 0;
}
