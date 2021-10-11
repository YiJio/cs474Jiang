#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string>								// to_string
#include <cstring>								// strcat

#include "image.h"
#include "WriteImage.h"
#include "ReadImage.h"
#include "ReadImageHeader.h"

int computeAverage(ImageType& image, int row, int col, int masksize, int max);

void getAverage(char fname[], ImageType& image, int masksize);

int computeGaussian(ImageType& image, int row, int col, int masksize, int max);

void getGaussian(char fname[], ImageType& image, int masksize);
