#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string>										// to_string
#include <cstring>										// strcat
#include <cmath>										// pow

#include "image.h"
#include "WriteImage.h"
#include "ReadImage.h"
#include "ReadImageHeader.h"

int imageQuantization(char fname[], ImageType& image, int bpp, int mode);
