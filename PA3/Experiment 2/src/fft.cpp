#include "fft.h"

/**
 * This function computes for DFT using FFT using bit-reverse swapping first and then computing for
 * n-samples for first half of samples. This is slightly adjusted from Dr. Bebis's provided code.
 * @param: data array of std::complex<float> values, n integer length, isign integer for type of
 * transform, integer row to start with
 * @pre: original data array values
 * @post: new data array values after transform
 * @return: none
 */
void fft(std::complex<float> data[], int n, int isign, int r) {
	// start at 1 for bit-reverse swapping, first and last element itself
	int j = n / 2;	
	for(int i = 1; i < (n - 1); i++) {
		// swap with values ahead, no self-swapping/double swapping
		if(j > i) { std::swap(data[i * r], data[j * r]); }
		int m = n / 2;
		// find first 0 from left moving right and set to 1
		while(m > 0 && j & m) {
			// XOR value if nonzero and shift right continue
			j ^= m;
			m >>= 1;
		}
		// OR to correct bit
		j |= m;
	}
	
	// fft compute M fourier transforms to n-sample fourier transform
	for(int M = 2; M <= n; M *= 2) {
		// inverse or forward depend on isign, start at (1+0i)
		double theta = isign * (2 * M_PI / M);
		std::complex<double> W_M = {cos(theta), sin(theta)};
		std::complex<double> W_Mu = {1, 0};
		// calculate first half of samples
		for(int u = 0; u < M / 2; u++) {
			for(int i = u; i < n; i += M) {
				int j = i + M / 2;
				std::complex<double> curr = std::complex<double>(data[j * r]) * W_Mu;
				// F(u+M/2)=F_even(u)-F_odd(u)*M_M^u
				data[j * r] = std::complex<double>(data[i * r]) - curr;
				// F_odd(u)*W_M^u
				data[i * r] += curr;
			}
			// W_M^{u+1}=W_M^u*W_M
			W_Mu *= W_M;
		}
	}
}

/**
 * This function computes for DFT using 2D FFT by the separability property. This works by going through
 * each column and row on the function and computing FFT for it. Finally, the data is adjusted/corrected
 * by dividing it from the square root of the total number of values after data has been computed for
 * each section.
 * @param: data array of std::complex<float> values, N, M integer sizes of array, isign integer for
 * type of transform
 * @pre: original data array values
 * @post: new data array values after transform
 * @return: none
 */
void fft2D(std::complex<float> data[], int N, int M, int isign) {
	// fft on columns
	for(int i = 0; i < N; i++) { fft(data + i * M, M, isign); }
	// fft on rows
	for(int i = 0; i < M; i++) { fft(data + i, N, isign, M); }
	// correct fft on factor
	for(int i = 0; i < N * M; i++) { data[i] *= 1 / sqrt(N * M); }
}

/**
 * This function generates a black image with a square at the center. It takes in arguments that tell
 * how big the square should be. The image size is set to 512x512.
 * @param: n, m integer sizes for square, Q integer for maximum quantization level
 * @pre: all original values in variables
 * @post: image written out
 * @return: none
 */
void generateImage(int n, int m, int Q) {
	// variables
	int N = 512, M = 512;
	ImageType image(N, M, Q);
	
	// fill image black with nxn center white
	for(int i = 0; i < N; i++) {
		for(int j = 0; j < M; j++) {
			if((i>(N/2-n/2)) && (i<=(N-N/2+n/2)) && (j>(M/2-m/2)) && (j<=(M-M/2+m/2))) {
				image.setPixelVal(i, j, 255);
			} else { image.setPixelVal(i, j, 0); }
		}
	}

	// file names and writing
	std::string fname = "../images/img" + std::to_string(n) + ".pgm";
	char *imageFile = new char[fname.length() + 1];
	strcpy(imageFile, fname.c_str());
	writeImage(imageFile, image);
	delete[] imageFile;
}

/**
 * This function takes in an image and basically just grabs the value from the image and then stores
 * it to the transform array to perform the 2D FFT. However, it also performs translate/shift on the
 * magnitude to the center of the frequency domain before the transforms happen, if specified. It then
 * passes the transformed values to the 2D FFT and also calls to compute for the image values.
 * @param: fname character array to write image, image ImageType reference, transform array of
 * std::complex<float> values, size integer, mode integer preference
 * @pre: original transform array values
 * @post: new transform array values after transform
 * @return: none
 */
void transformImage(ImageType& image, std::complex<float> transform[], int size, int mode) {
	// variables
	int M, N, Q, value;
	float curr;
	image.getImageInfo(N, M, Q);
	
	// perform shift if mode is 1
	for(int i = 0; i < N; i++) {
		for(int j = 0; j < M; j++) {
			if(mode == 1) {
				// translate/shift magnitude to center of frequency domain
				if((i + j) % 2 == 0) { curr = 1; }
				else { curr = -1; }
				image.getPixelVal(i, j, value);
				transform[i * M + j] = {(float)value * curr, 0};
			} else {
				image.getPixelVal(i, j, value);
				transform[i * M + j] = {(float)value, 0};
			}
		}
	}
	fft2D(transform, N, M, -1);
	
	// call computeImage to shift or no shift the image
	computeImage(transform, N, M, size, mode, false);
	computeImage(transform, N, M, size, mode, true);
}

/**
 * This function computes for the image values by first grabbing the magnitude values from the transform
 * array. A copy of the transform is made so that the original values do not get altered when trying to
 * use the same transform properties on shift and no shift from the main() function. The magnitude value
 * is also altered if log transform is specified. The value is then stored in a temporary data array to
 * find out the min and max of the array to help with normalization. After normalization, the data is
 * then set to the image and written out.
 * @param: transform array of std::complex<float> values, N, M integer sizes of image, size integer
 * for square, mode integer preference, l bool log preference
 * @pre: original values in variables
 * @post: transformed image written out
 * @return: none
 */
void computeImage(std::complex<float> transform[], int N, int M, int size, int mode, bool l) {
	// variables
	int Q = 255;
	float data[N][M], value, curr, min, max;
	std::complex<float> test[N * M];
	ImageType image(N, M, Q);
		
	// grab the magnitude values from test transform to data array
	std::copy(transform, transform + (N * M), test);
	for(int i = 0; i < N; i++) {
		for(int j = 0; j < M; j++) {
			value = std::abs(test[i * M + j]);
			// c*log(1+F(u,v))
			if(l) { value = 1 * log(1 + value); }
			data[i][j] = value;
		}
	}
	
	// find the min and max values of the current data array to normalize values later
	min = data[0][0];
	max = min;
	for(int i = 0; i < N; i++) {
		for(int j = 0; j < M; j++) {
			curr = data[i][j];
			min = std::min(min, curr);
			max = std::max(max, curr);
		}
	}
	
	// must be computed after above step! so that the max and min values are final
	for(int i = 0; i < N; i++) {
		for(int j = 0; j < M; j++) {
			curr = data[i][j];
			int newvalue = 255l * (curr - min) / (max - min);
			image.setPixelVal(i, j, newvalue);
		}
	}

	// file names and writing
	std::string c = "";
	if(mode == 0) { c = "_n"; }
	else { c = "_y"; }
	std::string v = "";
	if(!l) { v = "_shift"; }
	else { v = "_shiftLog"; }
	std::string fname = "../images/img" + std::to_string(size) + c + v + ".pgm";
	char *imageFile = new char[fname.length() + 1];
	strcpy(imageFile, fname.c_str());
	writeImage(imageFile, image);
	delete[] imageFile;
}
