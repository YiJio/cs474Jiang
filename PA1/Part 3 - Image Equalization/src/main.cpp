#include <iostream>
#include <fstream>
#include <cstdlib>

#include "image.h"
#include "HistogramEqualization.h"

int main(int argc, char *argv[]) {

	std::cout << "Main called..." << std::endl;

	int N = 256;
	int M = 256;
	int Q = 255;
	int L = 256;

	// instantiate image objects
	ImageType boat(N, M, Q);
	ImageType f16(N, M, Q);

	// probabilities of each image
	double boat_pr[L] = {0.0};
	double f16_pr[L] = {0.0};
	
	// sets probabilities and equalizes boat
	getHistogram("boat", boat, boat_pr);
	equalizeImage("boat", boat, boat_pr);

	// sets probabilities and equalizes f_16
	getHistogram("f_16", f16, f16_pr);
	equalizeImage("f_16", f16, f16_pr);

	return 0;
}
