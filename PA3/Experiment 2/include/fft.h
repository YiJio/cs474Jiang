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

void fft(std::complex<float> data[], int n, int isign, int r = 1);

void fft2D(std::complex<float> data[], int N, int M, int isign);

void generateImage(int n, int m, int Q);

void transformImage(ImageType& image, std::complex<float> transform[], int size, int mode);

void computeImage(std::complex<float> transform[], int N, int M, int size, int mode, bool l);
