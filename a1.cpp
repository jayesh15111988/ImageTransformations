#include <SImage.h>
#include <SImageIO.h>
#include <cmath>
#include <ctime>
#include <math.h>
#include <unistd.h>
#include <getopt.h>
#define PI 3.14
#define E 2.718281828

//This is the name of default file. If no input is given, this file name will be considered
string input_filename = "4.png";

double sigma_value = 1.5;
int factor_value = 4;
int m_value = 7;
int n_value = 7;
double thresh_value = 5;
int GAUSS_KERNEL_SIZE = 3;
int hough_radius = 25;
int hough_vote = 900;
using namespace std;
SDoublePlane gaussian_filter(const SDoublePlane &input, double sigma);
// The simple image class is called SDoublePlane, with each pixel represented as
// a double (floating point) type. This means that an SDoublePlane can represent
// values outside the range 0-255, and thus can represent squared gradient magnitudes,
// harris corner scores, etc.
//
// The SImageIO class supports reading and writing PNG files. It will read in
// a color PNG file, convert it to grayscale, and then return it to you in
// an SDoublePlane. The values in this SDoublePlane will be in the range [0,255].
//
// To write out an image, call write_png_file(). It takes three separate planes,
// one for each primary color (red, green, blue). To write a grayscale image,
// just pass the same SDoublePlane for all 3 planes. In order to get sensible
// results, the values in the SDoublePlane should be in the range [0,255].
//
// Below are two helper functions that overlay rectangles / circles
// on an image plane for visualization purpose.
// Draws a rectangle on an image plane, using the specified gray level value and line width.
//
void overlay_rectangle(SDoublePlane &input, int _top, int _left, int _bottom, int _right, double graylevel, int width)
{
    for(int w=-width/2; w<=width/2; w++) {
        int top = _top+w, left = _left+w, right=_right+w, bottom=_bottom+w;
        // if any of the coordinates are out-of-bounds, truncate them
        top = min( max( top, 0 ), input.rows()-1);
        bottom = min( max( bottom, 0 ), input.rows()-1);
        left = min( max( left, 0 ), input.cols()-1);
        right = min( max( right, 0 ), input.cols()-1);
        // draw top and bottom lines
        for(int j=left; j<=right; j++)
            input[top][j] = input[bottom][j] = graylevel;
        // draw left and right lines
        for(int i=top; i<=bottom; i++)
            input[i][left] = input[i][right] = graylevel;
    }
}
// Draws a circle on an image plane, given circle center coordinate and radius.
//
void overlay_circle(SDoublePlane &input, int x_center, int y_center, int radius, int color)
{
    int r2 = radius * radius;
    for (int x = -radius; x <= radius; x++) {
        int y = (int) (sqrt(r2 - x*x)+0.5);
        if(x_center+x >= 0 && x_center+x < input.cols()) {
            if(y_center+y >=0 && y_center+y < input.rows())
                input[y_center+y][x_center+x] = color;
            if(y_center-y >=0 && y_center-y < input.rows())
                input[y_center-y][x_center+x] = color;
        }
    }
}
// The rest of these functions are incomplete. These are just suggestions to
// get you started -- feel free to add extra functions, change function
// parameters, etc.
// Given an image plane of size h x w, returns a subsampled image of size (h/k) x (w/k)
//
SDoublePlane subsample(const SDoublePlane &input, const int k)
{
    int rows = input.rows(), cols = input.cols(), i,j;
    SDoublePlane output(ceil(input.rows()/k)+1, ceil(input.cols()/k)+1);
    for(i=0; i<rows; i+=k) {
        for(j=0; j<cols; j+=k) {
            output[i/k][j/k] = input[i][j];
        }
    }
    return output;
}
// Convolve an image with a separable convolution kernel
//
//Following convolution for 1-D matrices was referred from 2-D convolution from source
//http://www.songho.ca/dsp/convolution/convolution.html#cpp_conv1d
SDoublePlane convolve_separable(const SDoublePlane &input, const SDoublePlane &row_filter, const SDoublePlane &col_filter)
{
    int rows=input.rows(), cols=input.cols(),i,j,k;
    SDoublePlane output(rows,cols);
    int kernel_size = row_filter.cols();
    int mid_row = floor(kernel_size/2);
    //create a temporary placeholder
    double ** temp_holder = new double *[rows];
    for(i=0; i<rows; i++)
        temp_holder[i] = new double[cols];
    //convolve the image with the 1-d vertical kernel
    for(i=0; i<rows; i++) {
        for(j=0; j<cols; j++) {
            temp_holder[i][j] = 0;
            for(k=0; k<kernel_size; k++) {
                int mask_row = kernel_size - 1 - k;
                int image_row = i + k - mid_row;
                if(image_row >= 0 && image_row < rows && j >= 0 && j < cols) {
                    temp_holder[i][j] += double(input[image_row][j] * col_filter[mask_row][0]);
                }
            }
        }
    }
    //convolve the image with the 1-d horizontal kernel
    for(i=0; i<rows; i++) {
        for(j=0; j<cols; j++) {
            output[i][j] = 0;
            for(k=0; k<kernel_size; k++) {
                int mask_col = kernel_size - 1 - k;
                int image_col = j + k - mid_row;
                if(image_col >= 0 && image_col < cols && i >= 0 && i < rows) {
                    output[i][j] += temp_holder[i][image_col] * row_filter[0][mask_col];
                }
            }
        }
    }
    return output;
}
// Apply a Gaussian of the specified sigma to an image, and return the result
//
SDoublePlane gaussian_filter(const SDoublePlane &input, double sigma)
{
    int kernel_size = GAUSS_KERNEL_SIZE, rows = input.rows(), cols = input.cols(),i,j,k,l;
    SDoublePlane output(rows,cols), row_filter(1,kernel_size), col_filter(kernel_size,1);
    double const1 = 1/(2 * PI * sigma * sigma);
    double const2 = -1/(2 * sigma * sigma);
    int mid_row = floor(kernel_size/2);
    int x = -mid_row;
    for(i=0; i<kernel_size; i++) {
        col_filter[i][0] = const1 * pow(E, (const2 * (pow(x,2) + pow(i+x,2))));
    }
    int y = 0;
    for(i=0; i<kernel_size; i++) {
        row_filter[0][i] = const1 * pow(E, (const2 * (pow(i+x,2) + pow(y,2))));
    }
    double value = row_filter[0][0];
    for(i=0; i<kernel_size; i++) {
        row_filter[0][i] = row_filter[0][i]/value;
    }
    output = convolve_separable(input,row_filter, col_filter);
    return output;
}
// Apply a mean filter to an image, returns the result
//
SDoublePlane mean_filter(const SDoublePlane &input, int m, int n)
{
    SDoublePlane output(input.rows(), input.cols());
    int rows = input.rows(), cols = input.cols(), i, j, k;
    int row_kernel_size = m;
    int col_kernel_size = n;
    //create a temporary array for convolving the image with the kernel
    double ** temp = new double * [rows];
    for(int i=0; i<rows; i++)
        temp[i] = new double[cols];
    double ** temp1 = new double * [rows];
    for(int i=0; i<rows; i++)
        temp1[i] = new double[cols];
    double ** op = new double * [rows];
    for(int i=0; i<rows; i++)
        op[i] = new double[cols];
    for(i=0; i<rows; i++) {
        double sum = 0;
        for(j=0; j<cols; j++) {
            temp[i][j] = 0;
            sum += input[i][j];
            temp[i][j] = sum;
        }
    }
    int row_mid = row_kernel_size/2;
    for(i=0; i<rows; i++) {
        for(j=0; j<cols; j++) {
            if(j <= row_mid) {
                op[i][j] = temp[i][j+row_mid]/row_kernel_size;
            } else if(j > (cols-1-row_mid)) {
                op[i][j] = (temp[i][cols-1] - temp[i][j-row_mid-1])/row_kernel_size;
            } else {
                op[i][j] = (temp[i][j+row_mid] - temp[i][j-row_mid-1])/row_kernel_size;
            }
        }
    }
    int c=0;
    for(i=0; i<cols; i++) {
        double sum = 0;
        for(j=0; j<rows; j++) {
            temp1[j][i] = 0;
            sum += op[j][i];
            temp1[j][i] = sum;
            c++;
        }
    }
    int col_mid = col_kernel_size/2;
    int cs=0,cm=0,ce=0;
    for(i=0; i<cols; i++) {
        for(j=0; j<rows; j++) {
            output[j][i] = 0;
            if(j <= col_mid) {
                output[j][i] = temp1[j+col_mid][i]/col_kernel_size;
                cs++;
            } else if(j > (rows-1-col_mid)) {
                ce++;
                output[j][i] = (temp1[rows-1][i] - temp1[j-col_mid-1][i])/col_kernel_size;
            } else {
                cm++;
                output[j][i] = (temp1[j+col_mid][i] - temp1[j-col_mid-1][i])/col_kernel_size;
            }
        }
    }
    return output;
}
// Apply a sobel operator to an image, returns the result
//
SDoublePlane sobel_gradient_filter(const SDoublePlane &input, bool _gx)
{
    int rows = input.rows(), cols = input.cols(),i,j,k;
    SDoublePlane sobel_row(1,3);
    SDoublePlane sobel_col(3,1);
    SDoublePlane output(rows, cols);
    //vertical or horizontal
    if (_gx==false) {
        sobel_col[0][0] = -1;
        sobel_col[1][0] = -2;
        sobel_col[2][0] = -1;
        sobel_row[0][0] = -1;
        sobel_row[0][1] = 0;
        sobel_row[0][2] = 1;
    } else {
        sobel_col[0][0] = -1;
        sobel_col[1][0] = 0;
        sobel_col[2][0] = 1;
        sobel_row[0][0] = -1;
        sobel_row[0][1] = -2;
        sobel_row[0][2] = -1;
    }
    //convolution function
    output = convolve_separable(input, sobel_row, sobel_col);
    return output;
}
// Apply an edge detector to an image, returns the binary edge map
//
SDoublePlane find_edges(const SDoublePlane &gradient_x, SDoublePlane &gradient_y, double thresh)
{
    int rows = gradient_x.rows(), cols = gradient_y.cols(), i, j, k;
    SDoublePlane output(rows, cols);
    for(i=0; i<rows; i++) {
        for(int j=0; j<cols; j++) {
            output[i][j]=sqrt(pow(gradient_x[i][j],2)+pow(gradient_y[i][j],2));
            if (output[i][j]>=thresh)
                output[i][j] = 0;
            else
                output[i][j] = 255;
        }
    }
    for(i=0; i<rows; i++) {
        for(j=0; j<cols; j++) {
            if(i<rows-1 && i>0 && j<cols-1 && j>0) {
                if(output[i][j]>output[i+1][j+1] && output[i][j]>=output[i-1][j-1])
                    output[i][j]=0;
                else
                    output[i][j]=255;
            }
        }
    }
    return output;
}
// Apply a circle detector to an image, returns the detection result
//
double ***matrix(const int range1,const int range2,SDoublePlane &input,int radiushough,double wtthresh)
{
    int a,b;
    double ***temp=0;
    temp = new double**[input.rows()];
    for (int i = 0; i < input.rows(); ++i) {
        temp[i] = new double*[input.cols()];
        for (int j = 0; j < input.cols(); ++j)
            temp[i][j] = new double[range1];
    }
    for(int i=0; i<input.rows(); i++) {
        for(int j=0; j<input.cols(); j++) {
            for(int k=0; k<range1; k++) {
                temp[i][j][k]=0;
            }
        }
    }
    int rr;
    for(int i=0; i<input.rows(); i++) {
        for(int j=0; j<input.cols(); j++) {
            if(input[i][j]==0) {
                for(int r=0; r<range1; r++) {
                    if(r+range2==radiushough) {
                        for(double theta=0; theta<=2*PI; theta+=0.00174) {
                            a=(round(i-((r+range2)*cos(theta))));
                            b=(round(j-((r+range2)*sin(theta))));
                            if(a<input.rows() &&a>=0 && b>=0 && b<input.cols()) {
                                temp[a][b][r]++;
                            }
                        }
                    }
                }
            }
        }
    }
    for(int i=0; i<input.rows(); i++) {
        for(int j=0; j<input.cols(); j++) {
            for(int r=0; r<range1; r++) {
                if(j<input.cols()-1 && i<input.rows()-1) {
                    if((abs(temp[i][j][r]-temp[i][j+1][r])<300) || (abs(temp[i][j][r]-temp[i+1][j][r])<300)||(abs(temp[i][j][r]-temp[i+1][j-1][r])<300)||(abs(temp[i][j][r]-temp[i][j+2][r])<300)||(abs(temp[i][j][r]-temp[i+2][j][r])<300)||(abs(temp[i][j][r]-temp[i+2][j][r-2])<300)||(abs(temp[i][j][r]-temp[i+1][j+1][r-1])<300)) {
                        temp[i][j][r]=0;
                    }
                }
            }
        }
    }
    return temp;
}
//To find the circles in the image
SDoublePlane hough_circles(SDoublePlane &input,int radiushough,double wtthresh)
{
    SDoublePlane output(input.rows(),input.cols());
    int x,y;
    int range1=40;
    int range2=10;
    int *row=0;
    int *column=0;
    int *radius=0;
    int r=0,rowmat=0,colmat=0;
    radius = new int[range1*input.rows()*input.cols()];
    row=new int[range1*input.rows()*input.cols()];
    column=new int[range1*input.rows()*input.cols()];
    for(int i=0; i<input.rows()*input.cols()*range1; i++) {
        row[i]=0;
        column[i]=0;
        radius[i]=0;
    }
    int a,b;
    double ***temp=0;
    int range=10;
    temp=matrix(range1,range2,input,radiushough,wtthresh);
    int var=0;
    for(int i=0; i<input.rows(); i++) {
        for(int j=0; j<input.cols(); j++) {
            for(int k=0; k<range1; k++) {
                if(temp[i][j][k]>wtthresh) {
                    row[rowmat++]=i;
                    column[colmat++]=j;
                    radius[r++]=k;
                }
            }
        }
    }
    int val=0;
    for(int i=0; i<input.rows()*input.cols()*range1; i++) {
        if(row[i]==0) {
            val=i;
            break;
        }
        for(double theta=0; theta<=2*PI; theta+=0.00174) {
            x=(round(row[i]+(radius[i]+range2)*cos(theta)));
            y=(round(column[i]+(radius[i]+range2)*sin(theta)));
            if (x<input.rows() && y<input.cols()&& x>=0 && y>=0) {
                output[x][y]=255;
            }
        }
    }
    for(int i=0; i<input.rows(); i++) {
        for(int j=0; j<input.cols(); j++) {
            delete [] temp[i][j];
        }
        delete [] temp[i];
    }
    delete temp;
    return output;
}
void displayOptions()
{
    cout<<"The following are the available options for the Program."<<endl;
    cout<<"1) image <image to process>"<<endl;
    cout<<"2) sigma <sigma value of the gaussian kernel>"<<endl;
    cout<<"3) factor <factor for subsampling>"<<endl;
    cout<<"5) thresh <threshold for edge detection>"<<endl;
    cout<<"6) gauss_kernel_size <size of the gauss kernel>"<<endl;
    cout<<"7) hough_radius <radius of the circle to be detected>"<<endl;
    cout<<"8) hough_vote <vote threshold of the pixels in hough space>"<<endl;
    cout<<"Please Enter the options in the following order in the command line. If you skip, default value will be taken"<<endl;
}
int main(int argc, char *argv[])
{
    if(argc < 2) {
        cerr << "usage: " << argv[0] << " input_image" << endl;
        displayOptions();
        return 1;
    }
    if(argc < 3) {
        displayOptions();
    }

    if(argc >=2)
        input_filename = argv[1];
    if(argc >= 3)
        sigma_value = atof(argv[2]);
    if(argc >= 4)
        factor_value = atoi(argv[3]);
    if(argc >= 5)
        thresh_value = atof(argv[4]);
    if(argc >= 6)
        GAUSS_KERNEL_SIZE = atoi(argv[5]);
    if(argc >= 7)
        hough_radius = atoi(argv[6]);
    if(argc >= 8)
        hough_vote = atoi(argv[7]);

    printf("Converting... Please be patient...\n\n");
    int start_s=clock();
    printf("Reading input file \n\n");
    SDoublePlane input_image= SImageIO::read_png_file(input_filename.c_str());
    //out of the gray scale image
    printf("Converting to gray scale image\n\n");
    SImageIO::write_png_file("gray_scale_image.png", input_image, input_image, input_image);
    // subsample the input image
    printf("Subsampling gray scale image\n\n");
    SDoublePlane subsampled = subsample(input_image,factor_value);
    SImageIO::write_png_file("subsampled_image.png", subsampled, subsampled, subsampled);
    //Call to Gaussian filter to smooth the image
    SDoublePlane gaussed = gaussian_filter(subsampled,sigma_value);
    printf("Applying Gaussian filter to image\n\n");
    SImageIO::write_png_file("gaussed_image.png", gaussed , gaussed , gaussed);
    //mean filter
    printf("Calculating mean value....\n\n");
    SDoublePlane mean_filtered = mean_filter(input_image, m_value, n_value);
    SImageIO::write_png_file("mean_image.png", mean_filtered, mean_filtered, mean_filtered);
    // compute gradient magnitude map of the input image
    SDoublePlane gradient_x = sobel_gradient_filter(gaussed, true);
    printf("Applying Sobel filter in x direction\n\n");
    SImageIO::write_png_file("gradient_x_image.png", gradient_x, gradient_x, gradient_x);
    SDoublePlane gradient_y = sobel_gradient_filter(gaussed, false);
    printf("Applying Sobel filter in y direction\n\n");
    SImageIO::write_png_file("gradient_y_image.png", gradient_y, gradient_y, gradient_y);
    // find edges in the input image
    SDoublePlane edges = find_edges(gradient_x, gradient_y, thresh_value);
    printf("Applying edge detection algorithm\n\n");
    SImageIO::write_png_file("edges_image.png", edges, edges, edges);
    // detect circles in the input image, and visualize the result
    SDoublePlane overlay_plane = hough_circles(edges, hough_radius, hough_vote);
    printf("Applying circle detection algorithm\n\n");
    SImageIO::write_png_file("detected_image.png",overlay_plane,overlay_plane,overlay_plane);

    printf("Image Transformation Successfully Completed\n\n");
    int stop_s=clock();
    printf("Total time taken for image operations is %f \n\n",(stop_s-start_s)/double(CLOCKS_PER_SEC));
}
