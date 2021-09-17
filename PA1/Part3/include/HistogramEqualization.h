#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string>								// to_string
#include <cstring>								// strcat
#include <cmath>								// round
#include <numeric>								// accumulate

#include "image.h"
#include "WriteImage.h"
#include "ReadImage.h"
#include "ReadImageHeader.h"

void printHistogram(char fname[], double prob[], int size, std::string type);

void getHistogram(char fname[], ImageType& image, double pr[]);

void equalizeImage(char fname[], ImageType& image, double pr[]);
