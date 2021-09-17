#include <iostream>
#include <fstream>
#include <cstdlib>

#include "image.h"

int readImageHeader(char fname[], int& N, int& M, int& Q, bool& type);
