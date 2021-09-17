#include "ImageSampling.h"

int imageSampling(char fname[], ImageType& image, int sampleFactor) {
	// variables
	int M, N, Q, M_sub, N_sub, Q_sub;
	int value, N_count = 0, M_count = 0;
	bool type;

	// original file names and reading
	std::string oldfname = std::string(fname) + ".pgm";
	char oldImageFile[oldfname.length() + 1];
	strcpy(oldImageFile, oldfname.c_str());
	readImageHeader(oldImageFile, N, M, Q, type);
	readImage(oldImageFile, image);

	// subsampling properties
	M_sub = M / sampleFactor;
	N_sub = N / sampleFactor;
	//Q_sub = Q / sampleFactor;		keep quantization level?
	ImageType newImage(N_sub, M_sub, Q);

	// perform subsampling
	for(int i = 0; i < N; i += sampleFactor) {
		M_count = 0;
		for(int j = 0; j < M; j += sampleFactor) {
			image.getPixelVal(i, j, value);
			newImage.setPixelVal(N_count, M_count, value);
			M_count++;
		}
		N_count++;
	}

	// subsampled file names and writing
	std::string newfname = std::string(fname) + "_ss" + std::to_string(sampleFactor) + ".pgm";
	char newImageFile[newfname.length() + 1];
	strcpy(newImageFile, newfname.c_str());
	writeImage(newImageFile, newImage);

	return 1;
}

int imageSizing(char fname[], ImageType& image, int sizeFactor) {
	// variables
	int M, N, Q, M_sub, N_sub, Q_sub;
	int value, N_count = 0, M_count = 0;
	bool type;

	// original file names and reading
	std::string oldfname = std::string(fname) + ".pgm";
	char oldImageFile[oldfname.length() + 1];
	strcpy(oldImageFile, oldfname.c_str());
	readImageHeader(oldImageFile, N_sub, M_sub, Q, type);
	readImage(oldImageFile, image);

	// resizing properties
	M = M_sub * sizeFactor;
	N = N_sub * sizeFactor;
	double factor = 1 / (double)sizeFactor;
	ImageType newImage(N, M, Q);

	// perform resizing
	for(int i = 0; i < N; i++) {
		N_count = i * factor;
		for(int j = 0; j < M; j++) {
			M_count = j * factor;
			image.getPixelVal(N_count, M_count, value);
			newImage.setPixelVal(i, j, value);
		}
	}

	// resized file names and writing
	std::string newfname = std::string(fname) + "_rs.pgm";
	char newImageFile[newfname.length() + 1];
	strcpy(newImageFile, newfname.c_str());
	writeImage(newImageFile, newImage);

	return 1;
}
