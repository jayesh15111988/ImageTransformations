#include <SImage.h>
#include <SImageIO.h>
#include <cmath>
#include <math.h>
#define PI 3.14
#define E 2.718281828

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
  //Call to Gaussian filter to smooth the image
  double sigma;
  cout<<"enter the sigma value"<<endl;
  cin>>sigma;
  SDoublePlane input1 = gaussian_filter(input,sigma);
  SImageIO::write_png_file("gaussed.png", input1, input1, input1);
  
  int rows = input1.rows(), cols = input1.cols(), i,j;
  
  SDoublePlane output(ceil(input.rows()/k)+1, ceil(input.cols()/k)+1); 

  for(i=0;i<rows;i+=k){
		for(j=0;j<cols;j+=k){
			output[i/k][j/k] = input1[i][j];
		}
  }

  return output;
}

// Convolve an image with a separable convolution kernel
//
SDoublePlane convolve_separable(const SDoublePlane &input, const SDoublePlane &row_filter, const SDoublePlane &col_filter)
{
  int rows=input.rows(), cols=input.cols(),i,j,k;
  SDoublePlane output(rows,cols);
  int kernel_size = row_filter.cols();
  int mid_row = floor(kernel_size/2);
  
  //create a temporary placeholder
  double ** temp_holder = new double *[rows];
  for(i=0;i<rows;i++)
	temp_holder[i] = new double[cols];

  //convolve the image with the 1-d vertical kernel
  for(i=0;i<rows;i++){
	for(j=0;j<cols;j++){
		temp_holder[i][j] = 0;
		for(k=0;k<kernel_size;k++){
			int mask_row = kernel_size - 1 - k;
			int image_row = i + k - mid_row;
			if(image_row >= 0 && image_row < rows && j >= 0 && j < cols){
				temp_holder[i][j] += double(input[image_row][j] * col_filter[mask_row][0]);
			}
		}
		
  	}
  }

  //convolve the image with the 1-d horizontal kernel
  for(i=0;i<rows;i++){
	for(j=0;j<cols;j++){
		output[i][j] = 0;
		for(k=0;k<kernel_size;k++){
			int mask_col = kernel_size - 1 - k;
			int image_col = j + k - mid_row;
			if(image_col >= 0 && image_col < cols && i >= 0 && i < rows){
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
  int kernel_size = 7, rows = input.rows(), cols = input.cols(),i,j,k,l;
  SDoublePlane output(rows,cols), row_filter(1,kernel_size), col_filter(kernel_size,1);


  double const1 = 1/(2 * PI * sigma * sigma);
  double const2 = -1/(2 * sigma * sigma);
  int mid_row = floor(kernel_size/2);
  
  
  
  int x = -mid_row;
  for(i=0;i<kernel_size;i++){
	col_filter[i][0] =  const1 * pow(E, (const2 * (pow(x,2) + pow(i+x,2))));
  }
  
  int y = 0;
  for(i=0;i<kernel_size;i++){
 	row_filter[0][i] =  const1 * pow(E, (const2 * (pow(i+x,2) + pow(y,2))));
  }

  double value = row_filter[0][0];
  for(i=0;i<kernel_size;i++){
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
  for(int i=0;i<rows;i++)
	temp[i] = new double[cols];

  double ** temp1 = new double * [rows];
  for(int i=0;i<rows;i++)
	temp1[i] = new double[cols];

  double ** op = new double * [rows];
  for(int i=0;i<rows;i++)
	op[i] = new double[cols];

  for(i=0;i<rows;i++){
	double sum = 0;
	for(j=0;j<cols;j++){
		temp[i][j] = 0;
		sum += input[i][j];
		temp[i][j] = sum;
	}
  }
//cout<<"jayesh";
  int row_mid = row_kernel_size/2;

  for(i=0;i<rows;i++)
{

	for(j=0;j<cols;j++)
{
		if(j < row_mid)
{
			op[i][j] = temp[i][j+row_mid]/row_kernel_size;
		}else if(j > (cols-1-row_mid)){
			op[i][j] = (temp[i][cols-1] - temp[i][j-row_mid-1])/row_kernel_size;
		}else{
			op[i][j] = (temp[i][j+row_mid] - temp[i][j-row_mid-1])/row_kernel_size;
		}
	}
  }

//  cout<<temp[i-2][j-1]<<" "<<temp[i-2][j-2]<<" "<<temp[i-2][j-3]<<" "<<temp[i-2][j-4]<<" "<<temp[i-2][j-5]<<" "<<op[i-2][j-1]<<" "<<op[i-2][j-2]<<" "<<op[i-2][j-3]<<" "<<" "<<endl;
//cout<<"jayesh"<<endl;
  for(i=0;i<cols;i++){
	double sum = 0;
	for(j=0;j<rows;j++){
	temp1[j][i] = 0;
		sum += op[j][i];
		temp1[j][i] = sum;
	}
  }


  int col_mid = col_kernel_size/2;
  int cs=0,cm=0,ce=0;
  for(i=0;i<cols;i++){
	for(j=0;j<rows;j++){
		output[j][i] = 0;
		if(j < col_mid){
		//	output[j][i] = temp1[j+col_mid][i]/col_kernel_size;
//			cout<<"start"<<endl;
			cs++;
		}else if(j > (rows-1-col_mid)){
		//	output[j][i] = (temp1[rows-1][i] - temp1[j-col_mid-1][i])/col_kernel_size;
			cout<<"end"<<endl;
			ce++;
		}else{
cout<<j-col_mid-1<<endl;
			output[j][i] = (temp1[j+col_mid][i] - temp1[j-col_mid-1][i])/col_kernel_size;
			cout<<"mid"<<endl;
			cm++;
		}
	}
  }
  
  cout<<cs<<" "<<cm<<" "<<ce<<endl;
  cout<<output[j-1][i-1]<<" "<<output[j-2][i-1]<<" "<<output[j-3][i-1]<<" "<<output[j-4][i-1]<<" "<<output[j-5][i-1]<<" "<<temp1[j-1][i-1]<<" "<<temp1[j-2][i-1]<<" "<<temp1[j-3][i-1]<<" "<<" "<<endl;

  exit(1);
  // Implement a m*n 2D mean filter with 1-d filters
  
  return output;
}

// Apply a sobel operator to an image, returns the result
// 
SDoublePlane sobel_gradient_filter(const SDoublePlane &input, bool _gx)
{
  SDoublePlane output(input.rows(), input.cols());

  // Implement a sobel gradient estimation filter with 1-d filters
  

  return output;
}

// Apply an edge detector to an image, returns the binary edge map
// 
SDoublePlane find_edges(const SDoublePlane &input, double thresh=0)
{
  SDoublePlane output(input.rows(), input.cols());

  // Implement an edge detector of your choice,
  // use your sobel gradient operator to compute the gradient magnitude
  
  return output;
}

// Apply a circle detector to an image, returns the detection result
// 
SDoublePlane hough_circles(const SDoublePlane &input, const int x_bin, const int y_bin, const int r_bin)
{
  // Hough transform for circle detection
  // This skeleton code places a fixed cirlce at the image center. 
  
  SDoublePlane overlay_plane = input;
  overlay_circle(overlay_plane, input.cols()/2, input.rows()/2, 10, 255);
  return overlay_plane;
}



int main(int argc, char *argv[])
{
  if(!(argc == 2))
    {
      cerr << "usage: " << argv[0] << " input_image" << endl;
      return 1;
    }
  string input_filename(argv[1]);
  SDoublePlane input_image= SImageIO::read_png_file(input_filename.c_str());

  //SImageIO::write_png_file("gray_scale.png", input_image, input_image, input_image);
  
  /*// subsample the input image
  //cout<<"enter the factor by which you need to subsample"<<endl;
  //int factor;
  //cin>>factor;
  //SDoublePlane subsampled = subsample(input_image,factor);
  //SImageIO::write_png_file("subsampled.png", subsampled, subsampled, subsampled);*/

  // mean filter
 // cout<<"enter the size of the mXn kernel"<<endl;
  //int m,n;
  //cin>>m;
  //cin>>n;
   
     SDoublePlane mean_filtered = mean_filter(input_image, 3, 3);
  
  // compute gradient magnitude map of the input image
 /* SDoublePlane gradient_x = sobel_gradient_filter(input_image, true);
 // SDoublePlane gradient_y = sobel_gradient_filter(input_image, false);
 /// SDoublePlane gradmag(input_image.rows(), input_image.cols());
 // SImageIO::write_png_file("gradmag.png", gradmag, gradmag, gradmag);
  
  // find edges in the input image
  //SDoublePlane edges = find_edges(input_image);
 // SImageIO::write_png_file("edges.png", edges, edges, edges);
  
  // detect circles in the input image, and visualize the result
  //SDoublePlane overlay_plane = hough_circles(input_image, input_image.cols()/2, input_image.rows()/2, 10);
  //SImageIO::write_png_file("detected.png",input_image,overlay_plane,input_image);*/
}
