#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string>								// to_string
#include <cstring>								// strcat
#include <cmath>
#include <complex>

#include "image.h"
#include "WriteImage.h"
#include "ReadImage.h"
#include "ReadImageHeader.h"

int computeGaussian(ImageType& image, int row, int col, int masksize, int max);

void getGaussian(char fname[], ImageType& image, int masksize);

void fft(std::complex<float> data[], int n, int isign, int r = 1);

void fft2D(std::complex<float> data[], int N, int M, int isign);

void transformImage(char fname[], ImageType& image, std::complex<float> transform[]);

void getImage(char fname[], std::complex<float> transform[], int N, int M, int mode, bool l);

void removeNoise(std::complex<float> transform[], int N, int M, bool r);
