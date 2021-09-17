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

int closest(int a, int b, int val);

int match(int arr[], int n, int val);

int indexOf(int arr[], int n, int val);

void printHistogram(char fname[], double prob[], int size, std::string type);

void getHistogram(char fname[], ImageType& image, double pr[]);

void specifyImage(char fname[], ImageType& image, double pr[], double pz_s[]);
