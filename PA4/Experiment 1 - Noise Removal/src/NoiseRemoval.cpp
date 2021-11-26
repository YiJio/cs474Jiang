#include "NoiseRemoval.h"

/**
 * This function computes for gaussian smooth filtering of the image. Based on the mask size being
 * passed in, it will create a window that is of that mask size, which helps to determine the offset
 * and bounds of the window. It will make a window and then take the accumulated sum from all the
 * window's values multiplied by the gaussian mask value. This value is determined by the mask size,
 * which uses either mask7 or mask15. The sum is then averaged by the sum of all weights in the mask
 * and returned.
 * @param: image ImageType reference, integer values for original row and col for current pixel,
 * integer for max bound
 * @return: integer average
 */
int computeGaussian(ImageType& image, int row, int col, int masksize, int max) {	
	// gaussian mask 7x7
	int mask7[7][7] = {
		{1,1,2,2,2,1,1},
		{1,2,2,4,2,2,1},
		{2,2,4,8,4,2,2},
		{2,4,8,16,8,4,2},
		{2,2,4,8,4,2,2},
		{1,2,2,4,2,2,1},
		{1,1,2,2,2,1,1}
	};
	// gaussian mask 15x15
	int mask15[15][15] = {
		{2,2,3,4,5,5,6,6,6,5,5,4,3,2,2},
		{2,3,4,5,7,7,8,8,8,7,7,5,4,3,2},
		{3,4,6,7,9,10,10,11,10,10,9,7,6,4,3},
		{4,5,7,9,10,12,13,13,13,12,10,9,7,5,4},
		{5,7,9,11,13,14,15,16,15,14,13,11,9,7,5},
		{5,7,10,12,14,16,17,18,17,16,14,12,10,7, 5},
		{6,8,10,13,15,17,19,19,19,17,15,13,10,8,6},
		{6,8,11,13,16,18,19,20,19,18,16,13,11,8,6},
		{6,8,10,13,15,17,19,19,19,17,15,13,10,8,6},
		{5,7,10,12,14,16,17,18,17,16,14,12,10,7,5},
		{5,7,9,11,13,14,15,16,15,14,13,11,9,7,5},
		{4,5,7,9,10,12,13,13,13,12,10,9,7,5,4},
		{3,4,6,7,9,10,10,11,10,10,9,7,6,4,3},
		{2,3,4,5,7,7,8,8,8,7,7,5,4,3,2},
		{2,2,3,4,5,5,6,6,6,5,5,4,3,2,2}
	};

	// variables
	int k = row;
	int l = col;
	int o = masksize / 2;
	int window[masksize][masksize];
	
	// get window for current pixel
	for(int i = 0; i < masksize; i++, k++) {
		l = col;
		for(int j = 0; j < masksize; j++, l++) {
			if((k-o>=0 && k-o<max) && (l-o>=0 && l-o<max)) {
				image.getPixelVal(k-o, l-o, window[i][j]);
			}
			else { window[i][j] = 0; }
		}
	}
	
	// grab sums of gaussian and return average of it by summask
	int sum = 0;
	int summask = 0;
	for(int i = 0; i < masksize; i++) {
		for(int j = 0; j < masksize; j++) {
			if(masksize == 7) {
				sum += window[i][j] * mask7[i][j];
				summask += mask7[i][j];
			}
			else if(masksize == 15) {
				sum += window[i][j] * mask15[i][j];
				summask += mask15[i][j];
			}			
		}
	}
	int avg = sum / summask;
	
	return avg;
}

/**
 * This function takes in the image's info and then calls computeGaussian to compute the gaussian filter
 * of the current pixel's window and writes the new image out.
 * @param: fname character array to write image, image ImageType reference, integer mask size
 * @pre: all original values in variables
 * @post: gaussian image written
 * @return: none
 */
void getGaussian(char fname[], ImageType& image, int masksize) {
	// variables
	int M, N, Q, value;
	image.getImageInfo(N, M, Q);
	ImageType newImage(N, M, Q);

	// perform gaussian
	for(int i = 0; i < N; i++) {
		for(int j = 0; j < M; j++) {
			value = computeGaussian(image, i, j, masksize, N);
			newImage.setPixelVal(i, j, value);
		}
	}

	// gaussian file names and writing
	std::string newfname = "../images/" + std::string(fname) + "_gaussian_mask" + std::to_string(masksize) + ".pgm";
	char *newImageFile = new char[newfname.length() + 1];
	strcpy(newImageFile, newfname.c_str());
	writeImage(newImageFile, newImage);
	delete[] newImageFile;
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
 * This function removes noise using the band-reject filter. It will compute for the fft of an image and
 * then apply the H(u,v) filter function that is then multiplied to the transform. The inverse fft is
 * then taken. It will then call to create the new image.
 * @param: fname character array to write image, image ImageType reference, w float for width of band,
 * d0 float for radius
 * @return: none
 */
void computeBand(char fname[], ImageType& image, int method, float w, float d0) {
	// variables
	int M, N, Q;
	image.getImageInfo(N, M, Q);
	ImageType newImage(N, M, Q);
	float u, v, duv, d01, d02;
	std::complex<float> H;
	std::complex<float>* transform = new std::complex<float>[N * M];

	// center and then fft
	transformImage(image, transform, 1);
	
	// perform noise removal using band-reject
	for(int i = 0; i < N; i++) {
		for(int j = 0; j < M; j++) { 
			u = i-(N/2);
			v = j-(M/2);
			duv = sqrt((u*u)+(v*v));
			// ideal 0, butterworth 1, gaussian 2
			if(method == 0) {
				d01 = d0-(w/2);
				d02 = d0+(w/2);
				if(d01 <= duv && duv <= d02) { H = 0; }
				else { H = 1; }
			} else if(method == 1) {
				H = 1/(1+pow((duv*w)/((duv*duv)-(d0*d0)),2*10));
			} else if(method == 2) {
				H = 1-exp(-pow((duv*duv)-(d0*d0)/(duv*w),2));
			}			
			// perform multiplication of filter
			transform[i*M+j] *= H;
		}
	}
	
	// inverse transform
	fft2D(transform, N, M, 1);
	
	// generate the new image
	std::string newfname = std::string(fname) + "_band";
	char *imageFile = new char[newfname.length() + 1];
	strcpy(imageFile, newfname.c_str());
	getImage(imageFile, transform, N, M, false, 0);
}

/**
 * This function removes noise using the notch-reject filter. It will compute for the fft of an image and
 * then apply the H(u,v) filter function that is then multiplied to the transform. The inverse fft is
 * then taken. It will then call to create the new image.
 * @param: fname character array to write image, image ImageType reference, w float for width of band,
 * d0 float for radius, uk, vk float for notch locations
 * @return: none
 */
void computeNotch(char fname[], ImageType& image, int method, float w, float d0, float uk, float vk) {
	// variables
	int M, N, Q;
	image.getImageInfo(N, M, Q);
	ImageType newImage(N, M, Q);
	float u, v, d1uv, d2uv, d3uv, d4uv;
	std::complex<float> H;
	std::complex<float>* transform = new std::complex<float>[N * M];

	// center and then fft
	transformImage(image, transform, 1);
	
	// perform noise removal using notch-reject
	for(int i = 0; i < N; i++) {
		for(int j = 0; j < M; j++) { 
			u = i-(N/2);
			v = j-(M/2);
			d1uv = sqrt(pow(u-uk,2)+pow(v-vk,2));
			d2uv = sqrt(pow(u-uk,2)+pow(v+vk,2));
			d3uv = sqrt(pow(u+uk,2)+pow(v-vk,2));
			d4uv = sqrt(pow(u+uk,2)+pow(v+vk,2));
			if(d1uv <= d0 || d2uv <= d0 || d3uv <= d0 || d4uv <= d0) { H = 0; }
			else { H = 1; }
			// perform multiplication of filter
			transform[i*M+j] *= H;
		}
	}
	
	// inverse transform
	fft2D(transform, N, M, 1);
	
	// generate the new image
	std::string newfname = std::string(fname) + "_notch";
	char *imageFile = new char[newfname.length() + 1];
	strcpy(imageFile, newfname.c_str());
	getImage(imageFile, transform, N, M, false, 0);
}
