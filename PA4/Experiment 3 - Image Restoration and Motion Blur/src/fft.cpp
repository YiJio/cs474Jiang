#include "fft.h"

/**
 * This function 
 * @param: 
 * @pre: 
 * @post: 
 * @return: 
 */
float box_muller(float m, float s) {
	// mean m, standard deviation s
	float x1 = 0, x2 = 0, y1 = 0, w;
	static float y2;
	static int use_last = 0;

	if(use_last) {
		y1 = y2;
		use_last = 0;
	}
	else {
		while(true) {
			x1 = 2.0 * ((float)rand() / RAND_MAX) - 1.0;
			x2 = 2.0 * ((float)rand() / RAND_MAX) - 1.0;
			w = x1 * x1 + x2 * x2;
			if(w <= 1.0) { break; }
		}
		w = sqrt( (-2.0 * log( w ) ) / w );
		y1 = x1 * w;
		y2 = x2 * w;
		use_last = 1;
	}
	return (m + y1 * s);
}

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
 * it to the transform array to perform the 2D FFT. It also calls to compute for the image values.
 * @param: fname character array to write image, image ImageType reference, transform array of
 * std::complex<float> values
 * @pre: original transform array values
 * @post: new transform array values after transform
 * @return: none
 */
void transformImage(ImageType& image, std::complex<float> transform[], int mode) {
	// variables
	int M, N, Q, value;
	float curr;
	image.getImageInfo(N, M, Q);
	
	// perform the shifts if necessary
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
}

/**
 * This function computes for the image values by first grabbing the real values from the transform
 * array. A copy of the transform is made so that the original values do not get altered when trying to
 * use the same transform properties from the main() function. The real value is also altered if log
 * transform is specified. The value is then stored in a temporary data array to find out the min and
 * max of the array to help with normalization. After normalization, the data is then set to the image
 * and written out.
 * @param: fname character array to write image, transform array of std::complex<float> values,
 * N, M integer sizes of image, mode integer preference
 * @pre: original values in variables
 * @post: transformed image written out
 * @return: none
 */
void getImage(char fname[], std::complex<float> transform[], int N, int M, bool l) {
	// variables
	int Q = 255;
	float data[N][M], value, curr, min, max;
	std::complex<float> test[N * M];
	ImageType image(N, M, Q);
		
	// grab the magnitude values from test transform to data array
	std::copy(transform, transform + (N * M), test);
	for(int i = 0; i < N; i++) {
		for(int j = 0; j < M; j++) {
			if((i + j) % 2 == 0) { curr = 1; }
			else { curr = -1; }
			value = test[i * M + j].real() * curr;
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

/**
 * This function 
 * @param: 
 * @pre: 
 * @post: 
 * @return: 
 */
void blurImage(char fname[], ImageType& image, int blur) {
	// variables
	int M, N, Q, val;
	image.getImageInfo(N, M, Q);
	float u, v, c, sinc;
	float data[N][M];
	std::complex<float> H, G;
	std::complex<float>* huv = new std::complex<float>[N * M];
	std::complex<float>* transform = new std::complex<float>[N * M];
	float a = 0.1, b = 0.1, k = 0.25;

	// center and then fft
	transformImage(image, transform, 1);
	
	// blur image and add noise
	for(int i = 0; i < N; i++) {
		for(int j = 0; j < M; j++) { 
			u = i-(N/2);
			v = j-(M/2);
			c = 3.14*((u*a)+(v*b));
			// do this since c could be 0 which we don't want indeterminate!
			if(c == 0) { sinc = 1; }
			else { sinc = sin(c)/c; }
			// replace (1/c)*sin(c) with sinc
			H = sinc*exp(std::complex<float>(0,-c));
			huv[i*M+j] = H;
			// add noise with box muller method
			G = box_muller(0.0, (float)blur);
			// apply multiplication
			transform[i*M+j] = (transform[i*M+j] * H) + G;
		}
	}
	
	// inverse transform
	fft2D(transform, M, N, 1);
	
	// generate the new image
	std::string newfname = std::string(fname) + "_blur" + std::to_string(blur);
	char *imageFile = new char[newfname.length() + 1];
	strcpy(imageFile, newfname.c_str());
	getImage(imageFile, transform, N, M, false);
}

/**
 * This function 
 * @param: 
 * @pre: 
 * @post: 
 * @return: 
 */
void unblurImage(char fname[], ImageType& image, int mode, float d0, float k) {
	// variables
	int M, N, Q, val;
	image.getImageInfo(N, M, Q);
	float u, v, duv, c, sinc;
	float data[N][M];
	std::complex<float> H, G, B;
	std::complex<float>* huv = new std::complex<float>[N * M];
	std::complex<float>* transform = new std::complex<float>[N * M];
	float a = 0.1, b = 0.1, mag = 10;
	
	// center and then compute fft
	transformImage(image, transform, 1);
	
	// unblur image using inverse 0 or wiener 1
	for(int i = 0; i < N; i++) {
		for(int j = 0; j < M; j++) { 
			u = i-(N/2);
			v = j-(M/2);
			float c = 3.14*(u*a+v*b);
			// do this since c could be 0 which we don't want indeterminate!
			if(c == 0) { sinc = 1; }
			else { sinc = sin(c)/c; }
			// replace (1/c)*sin(c) with sinc
			std::complex<float> H = sinc*exp(std::complex<float>(0,-c));
			huv[i * M + j] = H;
			duv = sqrt((u*u)+(v*v));
			// butterworth
			B = 1/(1+pow(duv/d0,2*5));
			// gaussian
			G = exp(-(duv*duv)/(2*d0*d0));
			// perform filtering
			if(mode == 1) {
				transform[i*M+j] *= (1.0f/H) * (std::norm(H)/(std::norm(H)+k));
			} else {
				transform[i*M+j] *= (1.0f/H) * B;
			}
		}
	}
	
	// compute inverse fft
	fft2D(transform, M, N, 1);
	
	// generate the new image
	std::string m = "_Inverse_r" + std::to_string(d0);
	if(mode == 1) { m = "_Wiener_K" + std::to_string(k); }
	std::string newfname = std::string(fname) + m;
	char *imageFile = new char[newfname.length() + 1];
	strcpy(imageFile, newfname.c_str());
	getImage(imageFile, transform, N, M, false);
}
