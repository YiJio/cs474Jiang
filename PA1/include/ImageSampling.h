#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string>											// to_string
#include <cstring>										// strcat

#include "image.h"
#include "WriteImage.h"
#include "ReadImage.h"
#include "ReadImageHeader.h"

int imageSampling(char fname[], ImageType& image, int sampleFactor);

int imageSizing(char fname[], ImageType& image, int sizeFactor);
