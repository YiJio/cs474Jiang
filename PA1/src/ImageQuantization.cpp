#include "ImageQuantization.h"

int imageQuantization(char fname[], ImageType& image, int bpp, int mode) {
	// variables
	int M, N, Q, L, Q_new;
	int value;
	bool type;

	// original file names and reading
	std::string oldfname = std::string(fname) + ".pgm";
	char oldImageFile[oldfname.length() + 1];
	strcpy(oldImageFile, oldfname.c_str());
	readImageHeader(oldImageFile, N, M, Q, type);
	readImage(oldImageFile, image);

	// quantization properties
	L = pow(2, bpp);
	if(mode == 0) {
		// Q=L-1 (quantization reduced to bpp level)
		Q_new = L - 1;
	} else {
		// Q=255 (quantization keeping highest level)
		Q_new = 255;
	}
	ImageType newImage(N, M, Q_new);

	// perform quantization
	for(int i = 0; i < N; i++) {
		for(int j = 0; j < M; j++) {
			int current = 0;
			image.getPixelVal(i, j, value);
			if(mode == 0) { current = value / pow(2, 8 - bpp); }
			else { current = value * pow(2, 8 - bpp); }
			//std::cout << current << " ";
			newImage.setPixelVal(i, j, current);
		}
		//std::cout << std::endl;
	}
	//std::cout << std::endl;

	// quantized file names and writing
	std::string newfname = std::string(fname);
	if(mode == 0) { newfname += "_L" + std::to_string(L) + ".pgm"; }
	else { newfname += "_Q" + std::to_string(Q_new) + ".pgm"; }
	char newImageFile[newfname.length() + 1];
	strcpy(newImageFile, newfname.c_str());
	writeImage(newImageFile, newImage);

	return 1;
}
