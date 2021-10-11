#include "Correlation.h"

int computeCorrelation(ImageType& image, ImageType& mask, int row, int col, int Nm, int Mm, int N, int M) {
	int k = row;
	int l = col;
	int value;
	int on = Nm / 2;
	int om = Mm / 2;
	int window[Nm][Mm];
	
	for(int i = 0; i < Nm; i++, k++) {
		l = col;
		for(int j = 0; j < Mm; j++, l++) {
			if((k-on>=0 && k-on<N) && (l-om>=0 && l-om<M)) {
				image.getPixelVal(k-on, l-om, window[i][j]);
			}
			else { window[i][j] = 0; }
		}
	}
	int sum = 0;
	for(int i = 0; i < Nm; i++) {
		for(int j = 0; j < Mm; j++) {
			mask.getPixelVal(i, j, value);
			sum += window[i][j] * value;
		}
	}
	//int avg = sum / (Nm*Mm);
	return sum;
}

void getCorrelation(char fname[], ImageType& image, ImageType& mask) {
	// variables
	int M, N, Q, Mm, Nm, value;
	image.getImageInfo(N, M, Q);
	mask.getImageInfo(Nm, Mm, Q);
	ImageType newImage(N, M, Q);

	// perform correlation
	for(int i = 0; i < N; i++) {
		for(int j = 0; j < M; j++) {
			value = computeCorrelation(image, mask, i, j, Nm, Mm, N, M);
			//std::cout << "hi";
			newImage.setPixelVal(i, j, value);
		}
	}

	// correlation file names and writing
	std::string newfname = "../images/" + std::string(fname) + "_correlation.pgm";
	char newImageFile[newfname.length() + 1];
	strcpy(newImageFile, newfname.c_str());
	writeImage(newImageFile, newImage);
}
