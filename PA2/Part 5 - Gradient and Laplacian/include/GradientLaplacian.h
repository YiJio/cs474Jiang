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

int computePrewitt(ImageType& image, int row, int col, int max);

void getPrewitt(char fname[], ImageType& image);

int computeSobel(ImageType& image, int row, int col, int max);

void getSobel(char fname[], ImageType& image);

int computeLaplacian(ImageType& image, int row, int col, int max);

void getLaplacian(char fname[], ImageType& image);
