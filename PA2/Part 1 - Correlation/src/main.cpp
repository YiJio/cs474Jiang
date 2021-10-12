#include <iostream>
#include <fstream>
#include <cstdlib>

#include "image.h"
#include "Correlation.h"

int main(int argc, char *argv[]) {
	int N = 442, M = 288, Q = 255;
	int Nm = 83, Mm = 55;
	
	// read image and mask and perform correlation
	ImageType image(N, M, Q);
	ImageType mask(Nm, Mm, Q);
	readImage("../images/Image.pgm", image);
	readImage("../images/Pattern.pgm", mask);
	getCorrelation("Image", image, mask);

	return 0;
}
