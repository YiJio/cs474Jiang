#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string>								// to_string
#include <cstring>								// strcat

#include "image.h"
#include "WriteImage.h"
#include "ReadImage.h"
#include "ReadImageHeader.h"

void imageSampling(char fname[], ImageType& image, int sampleFactor);

void imageSizing(char fname[], ImageType& image, int sizeFactor);
