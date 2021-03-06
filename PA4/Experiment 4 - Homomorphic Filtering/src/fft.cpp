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
	// fft on rows
	for(int i = 0; i < N; i++) { fft(data + i * M, M, isign); }
	// fft on columns
	for(int i = 0; i < M; i++) { fft(data + i, N, isign, M); }
	// correct fft on factor
	for(int i = 0; i < N * M; i++) { data[i] *= 1 / sqrt(N * M); }
}

/**
 * This function takes in an image and basically just grabs the value from the image and then stores
 * it to the transform array to perform the 2D FFT. If the mode is specified to 1, then the transform
 * will shift to the center. It also calls to compute for the image values.
 * @param: fname character array to write image, image ImageType reference, transform array of
 * std::complex<float> values, mode integer preference
 * @pre: original transform array values
 * @post: new transform array values after transform
 * @return: none
 */
void transformImage(ImageType& image, std::complex<float> transform[], int mode) {
	// variables
	int M, N, Q, value;
	float shift;
	image.getImageInfo(N, M, Q);
	
	// perform the shifts if necessary
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
 * This function computes for the image values by first grabbing the real values from the transform
 * array. A copy of the transform is made so that the original values do not get altered when trying to
 * use the same transform properties from the main() function. Depending on the desired mode, the image
 * generated will either be the image or the spectrum. The value is also altered if log transform is
 * specified. The value is then stored in a temporary data array where the min and max of the array is
 * determined so that the values are normalized later. After normalization, the data is then set to the
 * image and written out.
 * @param: fname character array to write image, transform array of std::complex<float> values,
 * N, M integer sizes of image, l bool logarithmic preference, mode integer preference
 * @pre: original values in variables
 * @post: transformed image written out
 * @return: none
 */
void getImage(char fname[], std::complex<float> transform[], int N, int M, bool l, int mode) {
	// variables
	int Q = 255;
	float data[N][M];
	float value, shift, curr, min, max;
	std::complex<float>* test = new std::complex<float>[N * M];
	ImageType image(N, M, Q);
		
	// grab the magnitude values from test transform to data array
	std::copy(transform, transform + (N * M), test);
	for(int i = 0; i < N; i++) {
		for(int j = 0; j < M; j++) {
			shift = pow(-1, i+j);
			// image 0, spectrum 1
			if(mode == 0) { value = test[i*M+j].real() * shift; }
			else if(mode == 1) { value = std::abs(test[i*M+j]); }
			if(l) { value = 20 * log(1+value); }
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
			int newvalue = 255l*(curr-min)/(max-min);
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

/**
 * This function takes the image to perform homomorphic filtering. The components (illumance and
 * reflectance) have to be separated out first by using the log. These values are stored in the transform
 * array to be taken to compute for the fft. Next, the function will calculate the H(u,v) function and
 * apply multiplication with the transformed array. The inverse fft is then taken and the components are
 * then added back in by using the exponent. The values are also normalized to the range of [0,255] and
 * written out as an image.
 * @param: fname character array to write image, image ImageType reference, yH, yL float for gamma values
 * @return: none
 */
void homomorphicFilter(char fname[], ImageType& image, float yH, float yL) {
	// variables
	int M, N, Q, val;
	image.getImageInfo(N, M, Q);
	ImageType newImage(N, M, Q);
	float u, v, e, value, min, max, curr, shift;
	float data[N][M];
	std::complex<float> H;
	std::complex<float>* transform = new std::complex<float>[N * M];
	float c = 1, d0 = 1.8;
	
	// first need to log to separate illuminance and reflectance
	for(int i = 0; i < N; i++) {
		for(int j = 0; j < M; j++) {
			shift = ((i+j)%2 == 0) ? 1 : -1;
			image.getPixelVal(i, j, val);
			value = log((float)val);
			transform[i*M+j] = {value * shift, 0};
		}
	}
	
	// compute fft on log values
	fft2D(transform, N, M, -1);
	
	// then do the huv
	for(int i = 0; i < N; i++) {
		for(int j = 0; j < M; j++) {
			u = i-(N/2);
			v = j-(M/2);
			e = -c*((u*u+v*v)/(d0*d0));
			H = (yH-yL) * (1-exp(e)) + yL;
			transform[i*M+j] *= H;
		}
	}

	// compute inverse fft
	fft2D(transform, M, N, 1);
	
	// need to exponent back values
	for(int i = 0; i < M; i++) {
		for(int j = 0; j < N; j++) {
			shift = ((i+j)%2 == 0) ? 1 : -1;
			value = transform[i*N+j].real() * shift;
			data[i][j] = exp(value);
		}
	}
	
	// now do normalization by getting min and max first
	max = -1000000.0;
	min = 10000000.0;
	for(int i = 0; i < N; i++) {
		for (int j = 0; j < M; j++) {
			value = data[i][j];
			max = std::max(value, max);
			min = std::min(value, min);
		}
	}

	// set the new correct image values
	for(int i = 0; i < N; i++) {
		for(int j = 0; j < M; j++) {
			curr = data[i][j];
			int newvalue = 255*(double)(curr-min)/(double)(max-min);
			newImage.setPixelVal(i, j, newvalue);
		}
	}
	
	// generate the new image
	std::string mm = "";
	if(yH == 1.5 && yL == 0.5) { mm = "_H15_L05"; }
	else { mm = "_H" + std::to_string((int)yH) + "_L" + std::to_string((int)yL); }
	std::string newfname = "../images/" + std::string(fname) + mm + ".pgm";
	char *imageFile = new char[newfname.length() + 1];
	strcpy(imageFile, newfname.c_str());
	writeImage(imageFile, newImage);
	delete[] imageFile;
}
