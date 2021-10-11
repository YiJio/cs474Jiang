#include <iostream>
#include <fstream>
#include <cstdlib>

#include "image.h"
#include "HistogramSpecification.h"

int main(int argc, char *argv[]) {

	std::cout << "Main called..." << std::endl;

	int N = 256;
	int M = 256;
	int Q = 255;
	int L = 256;

	// instantiate image objects
	ImageType boat(N, M, Q);
	ImageType sf(N, M, Q);
	ImageType f16(N, M, Q);
	ImageType peppers(N, M, Q);

	// probabilities of each image
	double boat_pr[L] = {0.0};
	double sf_pr[L] = {0.0};
	double f16_pr[L] = {0.0};
	double peppers_pr[L] = {0.0};
	
	// sets probabilities and specifies boat
	getHistogram("boat", boat, boat_pr);
	getHistogram("sf", sf, sf_pr);
	specifyImage("boat", boat, boat_pr, sf_pr);

	// sets probabilities and specifies f_16
	getHistogram("f_16", f16, f16_pr);
	getHistogram("peppers", peppers, peppers_pr);
	specifyImage("f_16", f16, f16_pr, peppers_pr);

	return 0;
}
