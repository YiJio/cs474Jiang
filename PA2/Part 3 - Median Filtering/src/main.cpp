#include <iostream>
#include <fstream>
#include <cstdlib>

#include "image.h"
#include "MedianFiltering.h"

int main(int argc, char *argv[]) {
	int N = 256, M = 256, Q = 255;
	
	ImageType lenna(N, M, Q);
	readImage("../images/lenna.pgm", lenna);
	saltedPepper("lenna", lenna, 30);
	saltedPepper("lenna", lenna, 50);
	ImageType spLenna30(N, M, Q);
	ImageType spLenna50(N, M, Q);
	readImage("../images/lenna_X30.pgm", spLenna30);
	readImage("../images/lenna_X50.pgm", spLenna50);
	getMedian("lenna_X30", spLenna30, 7);
	getMedian("lenna_X30", spLenna30, 15);
	getMedian("lenna_X50", spLenna50, 7);
	getMedian("lenna_X50", spLenna50, 15);
	
	ImageType boat(N, M, Q);
	readImage("../images/boat.pgm", boat);
	saltedPepper("boat", boat, 30);
	saltedPepper("boat", boat, 50);
	ImageType spBoat30(N, M, Q);
	ImageType spBoat50(N, M, Q);
	readImage("../images/boat_X30.pgm", spBoat30);
	readImage("../images/boat_X50.pgm", spBoat50);
	getMedian("boat_X30", spBoat30, 7);
	getMedian("boat_X30", spBoat30, 15);
	getMedian("boat_X50", spBoat50, 7);
	getMedian("boat_X50", spBoat50, 15);
	
	return 0;
}
