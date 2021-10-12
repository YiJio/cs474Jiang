#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string>								// to_string
#include <cstring>								// strcat
#include <time.h>								// rand
#include <algorithm>							// sort

#include "image.h"
#include "WriteImage.h"
#include "ReadImage.h"
#include "ReadImageHeader.h"

int computeAverage(ImageType& image, int row, int col, int masksize, int max);

void getAverage(char fname[], ImageType& image, int masksize);

int computeMedian(ImageType& image, int row, int col, int masksize, int max);

void getMedian(char fname[], ImageType& image, int masksize);

void saltedPepper(char fname[], ImageType &image, int X);
