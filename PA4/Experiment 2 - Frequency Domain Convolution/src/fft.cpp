#include "fft.h"

// sobel fy mask
int sobel[3][3] = {
	{-1,0,1},
	{-2,0,2},
	{-1,0,1}
};

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
void transformImage(ImageType& image, std::complex<float> transform[], int mode) {
	// variables
	int M, N, Q, value;
	float shift;
	image.getImageInfo(N, M, Q);
	
	// perform shift if mode is 1
	for(int i = 0; i < N; i++) {
		for(int j = 0; j < M; j++) {
			if(mode == 1) {
				// translate/shift magnitude to center of frequency domain
				shift = pow(-1, i+j);
				image.getPixelVal(i, j, value);
				transform[i*M+j] = {(float)value * shift, 0};
			} else {
				image.getPixelVal(i, j, value);
				transform[i*M+j] = {(float)value, 0};
			}
		}
	}
	fft2D(transform, N, M, -1);
}

/**
 * This function computes the convolution in the spatial filter.
 * @param: fname character array to write image, image ImageType reference
 * @pre: original array values
 * @post: new array values and image written out
 * @return: none
 */
void spatialFilter(char fname[], ImageType& image) {
	// variables
	int M, N, Q, value, masksize = 3;
	float curr = 0, min, max;
	image.getImageInfo(N, M, Q);
	ImageType newImage(N, M, Q);
	float data[N * M];
	
	// in image, get window and mask
	for(int i = 0; i < N; i++) {
		for(int j = 0; j < M; j++) {
			float val = 0;
			for(int k = -1; k < 2; k++) {
				for(int l = -1; l < 2; l++) {
					// compute by mask and equation
					if(i+k<0 || i+k>=N || j+l<0 || j+l>=M) { val += 0; }
					else {
						image.getPixelVal(i+k, j+l, value);
						val += value * sobel[masksize/2 - k][masksize/2 - l];
					}
				}
			}
			// apply new value computed
			data[i*M+j] = val;
		}
	}
	
	// find min and max
	min = data[0];
	max = min;
    for(int i = 0; i < N; i++) {
        for(int j = 0; j < M; j++) {
            curr = data[i*M+j];
            min = std::min(min, curr);
            max = std::max(max, curr);
        }
    }
	
	// must be computed after above step! so that the max and min values are final
	for(int i = 0; i < N; i++) {
		for(int j = 0; j < M; j++) {
			curr = data[i * M + j];
			int newvalue = 255l*(curr-min)/(max-min);
			newImage.setPixelVal(i, j, newvalue);
		}
	}
	
	// image file names
	std::string newfname = "../images/" + std::string(fname) + "_spatial.pgm";
	char *newFile = new char[newfname.length() + 1];
	strcpy(newFile, newfname.c_str());
	writeImage(newFile, newImage);
	delete[] newFile;
}

/**
 * This function computes the convolution in the frequency filter.
 * @param: fname character array to write image, image ImageType reference, huv array of
 * std::complex<float> values
 * @pre: original transform array values
 * @post: new transform array values and image written out
 * @return: none
 */
void frequencyFilter(char fname[], ImageType& image, std::complex<float> huv[]) {
	// variables
	int M, N, Q, masksize = 3;
	float curr, value, max, min, shift;
	image.getImageInfo(N, M, Q);
	ImageType newImage(N, M, Q);
	float hxy[N][M], data[M * N];
	std::complex<float>* transform = new std::complex<float>[N * M];
	
	// compute fft on image (not centered)
	transformImage(image, transform, 0);

	// create h(x,y) and pad with zeroes
	for(int i = 0; i < N; i++) {
		for(int j = 0; j < M; j++) { hxy[i][j] = 0; }
	}

	// place sobel mask at center of h(x,y)
	for(int i = 0; i < masksize; i++) {
		for(int j = 0; j < masksize; j++) {
			int k = N/2-1 + i;
			int l = M/2-1 + j;
			hxy[k][l] = sobel[i][j];
		}
	}

	// compute fft of h(x,y)
	for(int i = 0; i < N; i++) {
		for(int j = 0; j < M; j++) {
			huv[i*M+j] = std::complex<float>(hxy[i][j], 0);
		}
	}
	fft2D(huv, N, M, -1);

	// apply complex multiplication and inverse fft to get g(x,y)
	for(int i = 0; i < N; i++) {
		for(int j = 0; j < M; j++) {
			shift = pow(-1, i + j);
			huv[i*M+j] = std::complex<float>(0, huv[i*M+j].imag()) * shift;
			transform[i*M+j] *= huv[i*M+j];
		}
	}
	fft2D(transform, N, M, 1);

	// set image values (not centered)
	for(int i = 0; i < N; i++) {
		for(int j = 0; j < M; j++) {
			data[i*M+j] = transform[i*M+j].real();
		}
	}

	// find min and max
	max = -1000000.0;
	min = 10000000.0;
	for(int i = 0; i < N; i++) {
		for (int j = 0; j < M; j++) {
			value = data[i*N+j];
			max = std::max(value, max);
			min = std::min(value, min);
		}
	}

	// must be computed after above step! so that the max and min values are final
	for(int i = 0; i < N; i++) {
		for(int j = 0; j < M; j++) {
			curr = data[i*N+j];
			int newvalue = 255*(double)(curr-min)/(double)(max-min);
			newImage.setPixelVal(i, j, newvalue);
		}
	}

	// image file names
	std::string newfname = "../images/" + std::string(fname) + "_frequency.pgm";
	char *newFile = new char[newfname.length() + 1];
	strcpy(newFile, newfname.c_str());
	writeImage(newFile, newImage);
	delete[] newFile;
}

/**
 * This function computes for the image spectrum values by first grabbing the real values from the
 * transform array. A copy of the transform is made so that the original values do not get altered when
 * trying to use the same transform properties from the main() function. The absolute values are taken
 * and altered if log transform is specified. The value is then stored in a temporary data array to find
 * out the min and max of the array to help with normalization. After normalization, the data is then
 * set to the image and written out.
 * @param: fname character array to write image, transform array of std::complex<float> values,
 * N, M integer sizes of image, l bool logarithmic preference
 * @pre: original values in variables
 * @post: new values in variables, transformed image written out
 * @return: none
 */
void getImageSpectrum(char fname[], std::complex<float> transform[], int N, int M, bool l) {
	// variables
	int Q = 255;
	float data[N][M], value, curr, min, max;
	std::complex<float>* test = new std::complex<float>[N * M];
	ImageType image(N, M, Q);
		
	// grab the magnitude values from test transform to data array
	std::copy(transform, transform + (N * M), test);
	for(int i = 0; i < N; i++) {
		for(int j = 0; j < M; j++) {
			value = std::abs(test[i * M + j]);
			if(l) { value = 20 * log(1 + value); }
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
	std::string newfname = "../images/" + std::string(fname) + ".pgm";
	char *imageFile = new char[newfname.length() + 1];
	strcpy(imageFile, newfname.c_str());
	writeImage(imageFile, image);
	delete[] imageFile;
}
