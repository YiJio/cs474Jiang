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

void transformImage(ImageType& image, std::complex<float> transform[], int mode);

void getImage(char fname[], std::complex<float> transform[], int N, int M, bool l, int mode);

void homomorphicFilter(char fname[], ImageType& image, float yH, float yL);
