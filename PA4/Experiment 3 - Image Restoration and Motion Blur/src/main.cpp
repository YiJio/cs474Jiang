#include <iostream>
#include <fstream>
#include <cstdlib>

#include "image.h"
#include "fft.h"

int main(int argc, char *argv[]) {
	int N = 256, M = 256, Q = 255;
	
	// generate lenna image to test
	ImageType lenna(N, M, Q);
	readImage("../images/lenna_muller100.pgm", lenna);
	std::complex<float>* test = new std::complex<float>[N * M];

	float a = 0.1, b = 0.1, k = 0.025, curr, value, min, max, u, v;
	int mode = 1;
	float data[N][M];
	ImageType newImage(N, M, Q);
	
	transformImage("lenna_muller100", lenna, test);
	getImage("lenna_muller100_as", test, N, M, true);
	
	for(int i = 0; i < N; i++) {
		for(int j = 0; j < M; j++) { 
			u = i-(N/2);
			v = j-(M/2);
			float blur = 3.14*((u*a)+(v*b) + 0.00001);
			//float h = (1 / (blur)) * sin(blur) * exp(-(sqrt((u*u)+(v*v))) * blur);			
			//float h = exp(-k*pow(pow(u,2)+pow(v,2),0.8333333333));
			//float h = exp(-k*pow((u*u) + (v*v), (float)5/6));
				//std::cout <<u << "," << v <<"," << h << " ";
			float c = 3.14*(u*a+v*b);
			float h = (1/c)*sin(c)*exp(-c);
			//std::cout << h << " ";
			
			float rad = 256;
			float mag = 10;
			float duv = sqrt((u*u) + (v*v));
			float B = 1 / (1 + pow(duv / rad, 2*mag));
			
			if(mode == 1) {
				test[i*M+j] *= h;
				//test[i*M+j] = (1/h) * (abs(h*h) / (abs(h*h) + k)) * test[i*M+j];
				//std::cout << test[i*M+j] << " ";
			} else {
				test[i*M+j] = std::abs(test[i*M+j])/h * B;
			}
		}
	}
	
	fft2D(test, M, N, 1);
	// set image values
	for(int i = 0; i < N; i++) {
		for(int j = 0; j < M; j++) {
			data[i][j] = test[i * M + j].real();
			std::cout << data[i][j] << " ";
		}
	}
	
	max = -1000000.0;
	min = 10000000.0;

	for(int i = 0; i < N; i++) {
		for (int j = 0; j < M; j++) {
			value = data[i][j];
			max = std::max(value, max);
			min = std::min(value, min);
		}
	}

	for(int i = 0; i < N; i++) {
		for(int j = 0; j < M; j++) {
			curr = data[i][j];
			int newvalue = 255 * (double)(curr - min) / (double)(max - min);
			newImage.setPixelVal(i, j, newvalue);
		}
	}
	
	std::string newfname = "../images/lenna_muller100_out.pgm";
	char *imageFile = new char[newfname.length() + 1];
	strcpy(imageFile, newfname.c_str());
	writeImage(imageFile, newImage);
	delete[] imageFile;
	
	
	return 0;
}
