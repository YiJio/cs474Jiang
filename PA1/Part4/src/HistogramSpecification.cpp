#include "HistogramSpecification.h"

/**
 * This function helps return which two values are closest to a target value.
 * @param: a integer for first value, b integer for second value, val integer for target value
 * @return: the integer value closest to val
 */
int closest(int a, int b, int val) {
	if(val - a >= b - val) { return b; }
	return a;
}

/**
 * This function finds the closest match to a value from an array.
 * @param: arr integer array for the range of values, n integer for size of array, val integer
 * for target value
 * @return: the integer value that is the best match in array
 */
int match(int arr[], int n, int val) {
	// check corner cases
	if(val <= arr[0]) { return arr[0]; }
	if(val >= arr[n - 1]) { return arr[n - 1]; }
	// binary search
	int i = 0, j = n, mid = 0;
	while(i < j) {
		mid = (i + j) / 2;
		if(val < arr[mid]) {
			if(mid > 0 && val > arr[mid - 1]) { return closest(arr[mid - 1], arr[mid], val); }
			j = mid;
		} else {
			if(mid < n - 1 && val < arr[mid + 1]) { return closest(arr[mid], arr[mid + 1], val); }
			i = mid + 1;
		}
	}
	return arr[mid];
}

/**
 * This function gets the index of a value from an array.
 * @param: arr integer array to check, n integer for size of array, val integer for target
 * value
 * @return: the integer value for the index
 */
int indexOf(int arr[], int n, int val) {
	int i = 0;
	while(i < n) {
		if(arr[i] == val) break;
		i++;
	}
	return i;
}

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
 * This function first equalizes a source image by computing s=T(r), then computes v=G(z),
 * and finally computes z=G^-1(s) with specification histogram. It also prints histogram
 * values in comma-delimited format and writes another image after everything.
 * @param: fname character array to write image, image ImageType reference,
 * existing pr double array to access the probabilities of each gray level, specified pz_s
 * double array to access probabilities of specified image
 * @return: none
 */
void specifyImage(char fname[], ImageType& image, double pr[], double pz_s[]) {
	// image properties
	int N, M, Q, value;
	image.getImageInfo(N, M, Q);
	ImageType specifiedImage(N, M, Q);

	// histogram properties
	int L = Q + 1;
	int s[L] = {0};					// s (array for s=T(r))
	int v[L] = {0};					// v (array for v=G(z))
	int z[L] = {0};					// z (array for z=G^-1(s))
	double ps[L] = {0.0};				// probabilities of s
	double pz_a[L] = {0.0};				// probabilities of z actual

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

	// compute v=G(z)
	for(int i = 0; i < L; i++) {
		double vgz = 0.0;
		if(i == 0) { vgz = pz_s[i] * Q; }
		else {
			double temp = 0.0;
			temp = std::accumulate(pz_s, pz_s + i + 1, temp);
			vgz = temp * Q;
		}
		int v_k = (int)round(vgz);
		v[i] = v_k;
	}

	// compute z=G^-1(s)
	for(int i = 0; i < L; i++) {
		int s_k = s[i];
		// get closest match of s_k in v to v_k
		int v_k = match(v, L, s_k);
		// grab index of v_k in z to z_k
		int z_k = indexOf(v, L, v_k);
		z[i] = z_k;
		// prob_r to prob_z_actual based on z_k
		if(pz_a[z_k] == 0) { pz_a[z_k] = pr[i]; }
		else { pz_a[z_k] += pr[i]; }
	}

	// print values to file to generate histogram with excel
	printHistogram(fname, pr, L, "r");
	printHistogram(fname, pz_s, L, "z_s");
	printHistogram(fname, pz_a, L, "z_a");

	// compute specification on images
	for(int i = 0; i < N; i++) {
		for(int j = 0; j < M; j++) {
			image.getPixelVal(i, j, value);
			specifiedImage.setPixelVal(i, j, z[value]);
		}
	}

	// equalized file names and writing
	std::string newfname = "../images/" + std::string(fname) + "_specified.pgm";
	char newImageFile[newfname.length() + 1];
	strcpy(newImageFile, newfname.c_str());
	writeImage(newImageFile, specifiedImage);
}
