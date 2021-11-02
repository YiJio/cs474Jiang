#include <iostream>
#include <fstream>
#include <iostream>
#include <cstdlib>
#include <algorithm>
#include <complex>

#include "fft.h"

int main(int argc, char *argv[]) {
	std::ofstream out;

	// part a
	std::complex<float> f[4] = {2.0, 3.0, 4.0, 4.0};
	// compute forward fft on f
	fft(f, 4, -1);
	out.open("../files/parta.txt");
	for(int i = 0; i < 4; i++) {
		f[i] = f[i] / std::complex<float>(4, 0);
		out << i << "," << f[i].real() << "," << f[i].imag() << "," << std::abs(f[i]) << "," << std::arg(f[i]) << std::endl;
	}
	out.close();
	// compute inverse fft on f to verify
	fft(f, 4, 1);
	for (int i = 0; i < 4; i++) { std::cout << "F'(F(f)):\t" << f[i].real() << "\t" << f[i].imag() << "\t" << std::abs(f[i]) << "\t" << std::arg(f[i]) << std::endl; }

	// part b
	int N = 128;
	float step = 1 / (float)N;
	std::complex<float> coswave[N];
	for(int i = 0; i < N; i++) {
		// f(x) = cos(2*pi*u*x), u=8, N=128
		coswave[i] = std::cos(2 * M_PI * 8 * (float)(i * step));
	}
	// perform shift on magnitude to center
	float shift = 1;
	for(int i = 0; i < N; i++) {
		coswave[i] = coswave[i] * shift;
		shift *= -1;
	}
	// compute forward fft on cosine samples
	fft(coswave, N, -1);
	out.open("../files/partb.txt");
	for(int i = 0; i < N; i++) {
		coswave[i] = coswave[i] / std::complex<float>(N, 0);
		out << i << "," << coswave[i].real() << "," << coswave[i].imag() << "," << std::abs(coswave[i]) << "," << std::arg(coswave[i]) << std::endl;
	}
	out.close();
	
	// part c
	std::ifstream in("../files/Rect_128.dat");
	std::complex<float> rect[N];
	// read data from rect file
	std::string line;
	double curr = 0;
	for(int i = 0; i < N; i++) {
    	in >> curr;
    	rect[i] = curr * std::pow(-1, i);    	
    }
    // compute forward fft on rect samples
    fft(rect, N, -1);
    out.open("../files/partc.txt");
    for(int i = 0; i < N; i++) {
    	rect[i] = rect[i] / std::complex<float>(N, 0);
    	out << i << "," << rect[i].real() << "," << rect[i].imag() << "," << std::abs(rect[i]) << "," << std::arg(rect[i]) << std::endl;
    }
    out.close();

	return 0;
}
