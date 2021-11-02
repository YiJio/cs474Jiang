#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string>								// to_string
#include <cstring>								// strcat
#include <cmath>								// sqrt, pow

#include "image.h"
#include "WriteImage.h"
#include "ReadImage.h"
#include "ReadImageHeader.h"

void computeGradient(ImageType& image, int row, int col, int max, int mask, int& x, int& y, int& mag);

void getGradient(char fname[], ImageType& image, int mask);

int computeLaplacian(ImageType& image, int row, int col, int max);

void getLaplacian(char fname[], ImageType& image);
