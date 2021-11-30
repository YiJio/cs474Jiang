#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string>								// to_string
#include <cstring>								// strcat
#include <cmath>
#include <complex>
#include <random>
#include <math.h>

#include "image.h"
#include "WriteImage.h"
#include "ReadImage.h"
#include "ReadImageHeader.h"

float box_muller(float m, float s);

void fft(std::complex<float> data[], int n, int isign, int r = 1);

void fft2D(std::complex<float> data[], int N, int M, int isign);

void transformImage(ImageType& image, std::complex<float> transform[], int mode);

void getImage(char fname[], std::complex<float> transform[], int N, int M, bool l, int mode);

void blurImage(char fname[], int muller);

void unblurImage(char fname[], int mode, float d0, float k);
