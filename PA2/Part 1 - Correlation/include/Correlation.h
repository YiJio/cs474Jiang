#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string>								// to_string
#include <cstring>								// strcat

#include "image.h"
#include "WriteImage.h"
#include "ReadImage.h"
#include "ReadImageHeader.h"

int computeCorrelation(ImageType& image, ImageType& mask, int row, int col, int Nm, int Mm, int N, int M);

void getCorrelation(char fname[], ImageType& image, ImageType& mask);
