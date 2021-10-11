#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string>								// to_string
#include <cstring>								// strcat

#include "image.h"
#include "WriteImage.h"
#include "ReadImage.h"
#include "ReadImageHeader.h"

int computeGaussian(ImageType& image, int row, int col, int masksize, int max);

void getGaussian(char fname[], ImageType& image, int masksize);

void getSharp(char fname[], ImageType& image, ImageType& mask, int k);

void getSharpMask(char fname[], ImageType& image, ImageType& lpImage);
