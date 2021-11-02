#include "fft.h"

void fft(std::complex<float> data[], int n, int isign, int r) {
	// variables
	int j = n / 2;
	
	for(int i = 1; i < (n - 1); i++) {
		if(j > i) { std::swap(data[i * r], data[j * r]); }
		int m = n / 2;
		while(m > 0 && j & m) {
			j ^= m;
			m >>= 1;
		}
		j |= m;
	}
	
	for(int M = 2; M <= n; M *= 2) {
		double theta = isign * (2 * M_PI / M);
		std::complex<double> W_M = {cos(theta), sin(theta)};
		std::complex<double> W_Mu = {1, 0};
		for(int u = 0; u < M / 2; u++) {
			for(int i = u; i < n; i += M) {
				int j = i + M / 2;
				std::complex<double> curr = std::complex<double>(data[j * r]) * W_Mu;
				data[j * r] = std::complex<double>(data[i * r]) - curr;
				data[i * r] += curr;
			}
			W_Mu *= W_M;
		}
	}
}
