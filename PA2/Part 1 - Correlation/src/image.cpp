#include <cstdlib>

#include "image.h"

ImageType::ImageType()
{
 N = 0;
 M = 0;
 Q = 0;

 pixelValue = NULL;
}

ImageType::ImageType(int tmpN, int tmpM, int tmpQ)
{
 int i, j;

 N = tmpN;
 M = tmpM;
 Q = tmpQ;

 pixelValue = new int* [M];
 for(i=0; i<M; i++) {
   pixelValue[i] = new int[N];
   for(j=0; j<N; j++)
     pixelValue[i][j] = 0;
 }
}

ImageType::ImageType(ImageType& oldImage)
{
 int i, j;

 N = oldImage.N;
 M = oldImage.M;
 Q = oldImage.Q;

 pixelValue = new int* [M];
 for(i=0; i<M; i++) {
   pixelValue[i] = new int[N];
   for(j=0; j<N; j++)
     pixelValue[i][j] = oldImage.pixelValue[i][j];
 }
}

ImageType::~ImageType()
{
 int i;

 for(i=0; i<M; i++)
   delete [] pixelValue[i];
 delete [] pixelValue;
}


void ImageType::getImageInfo(int& rows, int& cols, int& levels)
{
 rows = M;
 cols = N;
 levels = Q;
} 

void ImageType::setImageInfo(int rows, int cols, int levels)
{
 M= rows;
 N= cols;
 Q= levels;
} 

void ImageType::setPixelVal(int i, int j, int val)
{
 pixelValue[i][j] = val;
}

void ImageType::getPixelVal(int i, int j, int& val)
{
 val = pixelValue[i][j];
}
