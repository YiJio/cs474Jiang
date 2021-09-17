#include "HistogramEqualization.h"

/**
 * This function uses an array of probabilities for each gray level and prints it out in
 * a comma-delimited format that can be accessed and plotted on Excel.
 * @param: fname character array to write text filename, existing pr double array to
 * access the probabilities of each gray level, size integer for pr array, type string
 * to determine what histogram type this is
 * @post: new file with histogram values
 * @return: none
 */
void printHistogram(char fname[], double prob[], int size, std::string type) {
	std::string filename = "../files/" + std::string(fname) + "_histo_" + type + ".txt";
	std::ofstream out;
	out.open(filename);
	if(out) {
		for(int i = 0; i < size; i++) {
			out << i << "," << prob[i] << std::endl;
		}
	} else {
		std::cerr << "Error: file could not be opened." << std::endl;
		exit(1);
	}
	out.close();
}

/**
 * This function reads the image and stores the data into the image reference, as well
 * as stores the probabilities of each gray level in the image.
 * @param: fname character array to read image, image ImageType reference,
 * pr double array to store the probabilities of gray levels
 * @pre: all original values in variables
 * @post: image has data, pr array has data
 * @return: none
 */
void getHistogram(char fname[], ImageType& image, double pr[]) {
	// variables
	int M, N, Q;
	int value;
	bool type;

	// original file names and reading
	std::string oldfname = "../images/" + std::string(fname) + ".pgm";
	char oldImageFile[oldfname.length() + 1];
	strcpy(oldImageFile, oldfname.c_str());
	readImageHeader(oldImageFile, N, M, Q, type);
	readImage(oldImageFile, image);

	// histogram properties
	int L = Q + 1;
	int sum = 0;
	int freq[L] = {0};			// frequencies

	// find frequencies of each gray level
	for(int i = 0; i < N; i++) {
		for(int j = 0; j < M; j++) {
			int current = 0;
			image.getPixelVal(i, j, current);
			freq[current]++;
		}
	}
	// find p_r(k) based on frequencies
	sum = std::accumulate(freq, freq + L, sum);
	for(int i = 0; i < L; i++) {
		pr[i] = freq[i] / (double)sum;
	}
}

/**
 * This function equalizes a source image, prints histogram values in comma-delimited
 * format, and writes another image after equalization.
 * @param: fname character array to write image, image ImageType reference,
 * existing pr double array to access the probabilities of each gray level
 * @return: none
 */
void equalizeImage(char fname[], ImageType& image, double pr[]) {
	// image properties
	int N, M, Q, value;
	image.getImageInfo(N, M, Q);
	ImageType equalizedImage(N, M, Q);

	// histogram properties
	int L = Q + 1;
	int s[L] = {0};					// s (array for equalized image)
	double ps[L] = {0.0};				// probabilities of s

	// compute s=T(r)
	for(int i = 0; i < L; i++) {
		double str = 0.0;
		if(i == 0) { str = pr[i] * Q; }
		else {
			double temp = 0.0;
			temp = std::accumulate(pr, pr + i + 1, temp);
			str = temp * Q;
		}
		int s_k = (int)round(str);
		s[i] = s_k;
		// prob_r to prob_s based on s_k
		if(ps[s_k] == 0) { ps[s_k] = pr[i]; }
		else { ps[s_k] += pr[i]; }
	}

	// print values to file to generate histogram with excel
	printHistogram(fname, pr, L, "r");
	printHistogram(fname, ps, L, "s");

	// compute equalization on image
	for(int i = 0; i < N; i++) {
		for(int j = 0; j < M; j++) {
			image.getPixelVal(i, j, value);
			equalizedImage.setPixelVal(i, j, s[value]);
		}
	}

	// equalized file names and writing
	std::string newfname = "../images/" + std::string(fname) + "_equalized.pgm";
	char newImageFile[newfname.length() + 1];
	strcpy(newImageFile, newfname.c_str());
	writeImage(newImageFile, equalizedImage);
}
