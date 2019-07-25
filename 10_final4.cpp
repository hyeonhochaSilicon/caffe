#include <opencv2/opencv.hpp>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

using namespace cv;
using namespace std;

void menu(int *ftype, int *stride);

int*** SepImageToCh(int ***, Mat, int);
int*** PaddingImage(int ***, int ***, Mat, int, int);

int*** ConvolutionEdge(int ***pchRGB, Mat mimage, int stride, int channel, int padding);
int*** ConvolutionBlur3x3(int ***pchRGB, Mat mimage, int ftype, int stride, int channel, int padding);
int*** ConvolutionBlur5x5(int ***pchRGB, Mat mimage, int ftype, int stride, int channel, int padding);
int*** ConvolutionBlur5x5(int ***pchRGB, Mat mimage, int ftype, int stride, int channel, int padding);

int*** MaxPooling(int***chRGB, Mat mimage, int filterW, int channel, int stride);

double getActivation(const double sigma);

int main(void)
{   
    int i, j, k;
    int ***pchRGB;
    int ***conv_RGB;  

    Mat mimage;
    mimage = imread("test3.jpg", IMREAD_COLOR);
 
    int channel;
    printf("Input number of channels:");
    scanf("%d", &channel);
    
    int ftype = -1;
    int filterW = 0;
    int padding = 0;
    int stride = 0;

    menu(&ftype, &stride);

    int ***chRGB;
    chRGB = SepImageToCh(chRGB, mimage, channel);                                //common function
    if(ftype < 0 || ftype > 3)
    {
        printf("please re select.\n");
    }
    else if(ftype == 3)
    {
        printf("필터 크기를 입력해주세요:");
        scanf("%d", &filterW);
        conv_RGB = MaxPooling(chRGB, mimage, filterW, channel, stride);
    }
    else
    {
        printf("Input number of paddings:");
        scanf("%d", &padding);
        printf("\n\n");
           
        pchRGB = PaddingImage(pchRGB, chRGB, mimage, channel, padding);
        
        switch (ftype)
        {
            case 0:
                conv_RGB = ConvolutionEdge(pchRGB, mimage, stride, channel, padding);
                break;
            
            case 1:
                conv_RGB = ConvolutionBlur3x3(pchRGB, mimage, ftype, stride, channel, padding);
                break;
            
            case 2:
                conv_RGB = ConvolutionBlur5x5(pchRGB, mimage, ftype, stride, channel, padding);
                break;

            default:
                break;
        }
    }    
    
    //free
    clock_t begin = clock();
    for(i = 0; i < channel; i++)
    {
        for(j = 0; j < mimage.rows; j++)
        {
            free(*(*(chRGB+i)+j));     
        }
        free(*(chRGB+i));
    }
    free(chRGB);

    if(ftype == 3)
    {
       for(i = 0; i < channel; i++)
       {
          for(j = 0; j < (((mimage.rows - filterW) / stride) + 1); j++)
          {
            free(*(*(conv_RGB+i)+j));     
          }
       }
       free(conv_RGB);
    }
    else
    {
        for(i = 0; i < channel; i++)
        {
           for(j = 0; j < mimage.rows + padding*2; j++)
           {
              free(*(*(pchRGB+i)+j));     
           }
        }
        free(pchRGB);
        
        switch(ftype)
        {
            case 0:
                for(i = 0; i < channel; i++)
                {
                    for(j = 0; j < ((mimage.rows - 3 + 2*padding) / stride) + 1; j++)
                    {
                        free(*(*(conv_RGB+i)+j));     
                    }
                }
                free(conv_RGB);
                break;
                
            case 1:
                for(i = 0; i < channel; i++)
                {
                    for(j = 0; j < ((mimage.rows - 3 + 2*padding) / stride) + 1; j++)
                    {
                        free(*(*(conv_RGB+i)+j));     
                    }
                }
                free(conv_RGB);
                break;
                
            case 2:
                for(i = 0; i < channel; i++)
                {
                    for(j = 0; j < ((mimage.rows - 5 + 2*padding) / stride) + 1; j++)
                    {
                        free(*(*(conv_RGB+i)+j));     
                    }
                }
                free(conv_RGB);
                break;

            default:
                break;
        }
    }

    clock_t end = clock();
    double elapsed_secs = double(end - begin);
    printf("free time: %lf\n\n", elapsed_secs / CLOCKS_PER_SEC);

    return 0;
}

void menu(int *ftype, int *stride)
{
    printf("0.Edge detect\n");
    printf("1.blur3x3\n");
    printf("2.blur5x5\n");
    printf("3.max_pooling\n");
    
    printf("Select menu:");
    scanf("%d", &(*ftype));
    
    printf("\n");
    printf("Input stride:");
    scanf("%d", &(*stride));
}

int*** SepImageToCh(int ***chRGB, Mat mimage, int channel)
{
    clock_t begin = clock();

    int k, i, j;

    chRGB = (int ***)malloc(channel*sizeof(int **));
    for(i = 0; i < channel; i++)
    {
        *(chRGB + i) = (int **)malloc(mimage.rows * sizeof(int *));
        for(j = 0; j < mimage.rows; j++)
        {
           *(*(chRGB +i) + j) = (int *)malloc(mimage.cols * sizeof(int));
        }
    }

    for(k = 0; k < channel; k++)
    {
       for(i = 0; i < mimage.rows; i++)
       {
           for(j = 0; j < mimage.cols; j++)
           {
              chRGB[k][i][j] = mimage.at<cv::Vec3b>(i, j)[k];
           }
       }
    }

    clock_t end = clock();
    double elapsed_secs = double(end - begin);
    printf("pixel seperate time: %lf\n\n", elapsed_secs / CLOCKS_PER_SEC);

    return chRGB;
}

int*** PaddingImage(int ***pchRGB, int ***chRGB, Mat mimage, int channel, int padding)
{
    clock_t begin = clock();

    int k, i, j;
    
    int rows = (mimage.rows + padding*2);
    int cols = (mimage.cols + padding*2);

    pchRGB = (int ***)malloc(channel*sizeof(int **));
    for(i = 0; i < channel; i++)
    {
        *(pchRGB + i) = (int **)malloc( rows * sizeof(int *) );
        for(j = 0; j < rows; j++)
        {
           *(*(pchRGB +i) + j) = (int *)malloc( cols * sizeof(int) );
        }
    }

    for(k = 0; k < channel; k++)
    {
       for(i = 0; i < rows; i++)
       {
          for(j = 0; j < cols; j++)
          {
             if(i >= padding && i < mimage.rows + padding)
             {
                if(j >= padding  && j <= mimage.cols + padding - 1)
                {
                    pchRGB[k][i][j] = chRGB[k][i-padding][j-padding];
                }
             } 
          }
       }
    } 

    clock_t end = clock();
    double elapsed_secs = double(end - begin);
    printf("padding time: %lf\n\n", elapsed_secs / CLOCKS_PER_SEC );

    return pchRGB;
}

       
int*** ConvolutionEdge(int ***pchRGB, Mat mimage, int stride, int channel, int padding)
{
    clock_t begin = clock();

    int ***conv_RGB;
    int k, i, j;
    int row;
    int col;
    
    int edge[3][3] = { 
                       {-1, -1, -1},
                       {-1, 8, -1},
                       {-1, -1, -1}
                     };
    
    int conv_rows = ((mimage.rows - 3 + 2*padding) / stride) + 1;
    int conv_cols = ((mimage.cols - 3 + 2*padding) / stride) + 1;

    conv_RGB = (int ***)malloc(channel * sizeof(int **));
    for (k = 0; k < channel; k++)
    {
       *(conv_RGB + k) = (int **)malloc(conv_rows * sizeof(int *));
       for (col = 0; col < conv_rows; col++)
       {
          *(*(conv_RGB + k) + col) = (int *)malloc(conv_cols * sizeof(int));
       }
    }
    
    for (k = 0; k < channel; k++)
    {
        for (row = 0; row < conv_rows; row++)
        {
            for (col = 0; col < conv_cols; col++)
	    {
	        for (i = 0; i < 3; i++)
	        {
	            for (j = 0; j < 3; j++)
                    {
                        conv_RGB[k][row][col] += pchRGB[k][i + stride * row][j + stride * col] * edge[i][j];
                    }
                }
            }
        }
    }    
   

    for(k = 0; k < channel; k++)
    {
       for(row = 0; row < conv_rows; row++)
       {
          for(col = 0; col < conv_cols; col++)
          {
             if (conv_RGB[k][row][col] < 0)
             {
	        conv_RGB[k][row][col] = 0;
	     }
	     else if (conv_RGB[k][row][col] > 255)
             {
	        conv_RGB[k][row][col] = 255;
	     }           
          }
       }
    }

    Mat picture(conv_rows, conv_cols, mimage.type());
    for (k = 0; k < channel; k++)
    {
       for (row = 0; row < conv_rows; row++)
       {
	  for (col = 0; col < conv_cols; col++)
	  {
	     picture.at<cv::Vec3b>(row, col)[k] = conv_RGB[k][row][col];
	  }
       }
    }
   
    imwrite("x.jpg", picture);

    clock_t end = clock();
    double elapsed_secs = double(end - begin);
    printf("Conv_Edge time: %lf\n\n", elapsed_secs / CLOCKS_PER_SEC);

    return conv_RGB;
}

int*** ConvolutionBlur3x3(int ***pchRGB, Mat mimage, int ftype, int stride, int channel, int padding)
{
    clock_t begin = clock();

    int ***conv_RGB;
    int k, i, j;
    int row;
    int col;
    int convsize;
    
    double blur3x3[3][3] = { 
                             {0.0625, 0.125, 0.0625},
                             {0.125, 0.25, 0.125},
                             {0.0625, 0.125, 0.0625}
                           };
    
    int conv_rows = ((mimage.rows - 3 + 2*padding) / stride) + 1;
    int conv_cols = ((mimage.cols - 3 + 2*padding) / stride) + 1;

    conv_RGB = (int ***)malloc(channel * sizeof(int **));
    for (i = 0; i < channel; i++)
    {
       *(conv_RGB + i) = (int **)malloc(conv_rows * sizeof(int *));
       for (j = 0; j < conv_rows; j++)
       {
          *(*(conv_RGB + i) + j) = (int *)malloc(conv_cols * sizeof(int));
       }
    }
    
    for (k = 0; k < channel; k++)
    {
        for (row = 0; row < conv_rows; row++)
        {
            for (col = 0; col < conv_cols; col++)
	    {
	        for (i = 0; i < 3; i++)
	        {
	            for (j = 0; j < 3; j++)
                    {
                        conv_RGB[k][row][col] += pchRGB[k][i + stride * row][j + stride * col] * blur3x3[i][j];
                    }
                }
            }
        }
    }    
   

    for(k = 0; k < channel; k++)
    {
       for(i = 0; i < conv_rows; i++)
       {
          for(j = 0; j < conv_cols; j++)
          {
             if (conv_RGB[k][i][j] < 0)
             {
	        conv_RGB[k][i][j] = 0;
	     }
	     else if (conv_RGB[k][i][j] > 255)
             {
	        conv_RGB[k][i][j] = 255;
	     }           
          }
       }
    }

    Mat picture(conv_rows, conv_cols, mimage.type());
    for (k = 0; k < channel; k++)
    {
       for (i = 0; i < conv_rows; i++)
       {
	  for (j = 0; j < conv_cols; j++)
	  {
	     picture.at<cv::Vec3b>(i, j)[k] = conv_RGB[k][i][j];
	  }
       }
    }
   
    imwrite("x.jpg", picture);

    clock_t end = clock();
    double elapsed_secs = double(end - begin);
    printf("Conv_Blur3x3: %lf\n\n", elapsed_secs / CLOCKS_PER_SEC);

    return conv_RGB;
}

int*** ConvolutionBlur5x5(int ***pchRGB, Mat mimage, int ftype, int stride, int channel, int padding)
{
    clock_t begin = clock();

    int ***conv_RGB;
    int k, i, j;
    int row;
    int col;
    int convsize;
    
    double blur5x5[5][5] = { 
                             {0.0039, 0.0156, 0.0234, 0.0156, 0.0039},
                             {0.0156, 0.0625, 0.0937, 0.0625, 0.0156},
                             {0.0234, 0.0937, 0.1406, 0.0937, 0.0234},
                             {0.0156, 0.0625, 0.0937, 0.0625, 0.0156},
                             {0.0039, 0.0156, 0.0234, 0.0156, 0.0039}                           
                           };
    
    int conv_rows = ((mimage.rows - 5 + 2*padding) / stride) + 1;
    int conv_cols = ((mimage.cols - 5 + 2*padding) / stride) + 1;

    conv_RGB = (int ***)malloc(channel * sizeof(int **));
    for (i = 0; i < channel; i++)
    {
       *(conv_RGB + i) = (int **)malloc(conv_rows * sizeof(int *));
       for (j = 0; j < conv_rows; j++)
       {
          *(*(conv_RGB + i) + j) = (int *)malloc(conv_cols * sizeof(int));
       }
    }
    
    for (k = 0; k < channel; k++)
    {
        for (row = 0; row < conv_rows; row++)
        {
            for (col = 0; col < conv_cols; col++)
	    {
	        for (i = 0; i < 5; i++)
	        {
	            for (j = 0; j < 5; j++)
                    {
                        conv_RGB[k][row][col] += pchRGB[k][i + stride * row][j + stride * col] * blur5x5[i][j];
                    }
                }
            }
        }
    }    
   

    for(k = 0; k < channel; k++)
    {
       for(i = 0; i < conv_rows; i++)
       {
          for(j = 0; j < conv_cols; j++)
          {
             if (conv_RGB[k][i][j] < 0)
             {
	        conv_RGB[k][i][j] = 0;
	     }
	     else if (conv_RGB[k][i][j] > 255)
             {
	        conv_RGB[k][i][j] = 255;
	     }           
          }
       }
    }

    Mat picture(conv_rows, conv_cols, mimage.type());
    for (k = 0; k < channel; k++)
    {
       for (i = 0; i < conv_rows; i++)
       {
	  for (j = 0; j < conv_cols; j++)
	  {
	     picture.at<cv::Vec3b>(i, j)[k] = conv_RGB[k][i][j];
	  }
       }
    }
   
    imwrite("x.jpg", picture);

    clock_t end = clock();
    double elapsed_secs = double(end - begin);
    printf("Conv_Blur5x5 time: %lf\n\n", elapsed_secs / CLOCKS_PER_SEC);

    return conv_RGB;
}

int*** MaxPooling(int ***chRGB, Mat mimage, int filterW, int channel, int stride)
{
   clock_t begin = clock();

   int ***conv_RGB;
   int i, j, k; 
   int row, col;
   int comp = 0;

   int pool_rows = (((mimage.rows - filterW) / stride) + 1);
   int pool_cols = (((mimage.cols - filterW) / stride) + 1);

    conv_RGB = (int ***)malloc(channel * sizeof(int **));
    for (i = 0; i < channel; i++)
    {
       *(conv_RGB + i) = (int **)malloc(pool_rows * sizeof(int *));
       for (j = 0; j < pool_rows; j++)
       {
          *(*(conv_RGB + i) + j) = (int *)malloc(pool_cols * sizeof(int));
       }
    }

    for (k = 0; k < channel; k++)
    {
       for (row = 0; row < pool_rows; row++)
       {
          for (col = 0; col < pool_cols; col++)
	  {
	     for (i = 0; i < filterW; i++)
	     {
	        for (j = 0; j < filterW; j++)
	        {
                   comp = chRGB[k][i + stride * row][j + stride * col];
                   if( comp > conv_RGB[k][row][col])
                   {
                      conv_RGB[k][row][col] = comp;
                   }
	        }
	     }
             
             comp = 0;
          }
       }
    }

    for(k = 0; k < channel; k++)
    {
       for(i = 0; i < pool_rows; i++)
       {
          for(j = 0; j < pool_cols; j++)
          {
             if (conv_RGB[k][i][j] < 0)
             {
	        conv_RGB[k][i][j] = 0;
	     }
	     else if (conv_RGB[k][i][j] > 255)
             {
	        conv_RGB[k][i][j] = 255;
	     }           
          }
       }
    }

    Mat picture(pool_rows, pool_cols, mimage.type());
    for (k = 0; k < channel; k++)
    {
       for (i = 0; i < pool_rows; i++)
       {
	  for (j = 0; j < pool_cols; j++)
	  {
	     picture.at<cv::Vec3b>(i, j)[k] = conv_RGB[k][i][j];
	  }
       }
    }
   
    imwrite("x.jpg", picture);

    clock_t end = clock();
    double elapsed_secs = double(end - begin);
    printf("Max_Pooling time: %lf\n\n", elapsed_secs / CLOCKS_PER_SEC);

    return conv_RGB;
}

double getActivation(const double sigma)
{
    return (1 / (1 + exp(-sigma)));
}

