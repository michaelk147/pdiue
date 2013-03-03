#include "opencv2/imgproc/imgproc.hpp"
#include "opencv2/highgui/highgui.hpp"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include "scalespacefuns.h"

using namespace cv;

/// Global variables
Mat src, dst, tmp;
char* window_name = "Pyramids Demo";


/**
 * @function main
 */
int main( int argc, char** argv )
{
  /// General instructions
  printf( "\n Zoom In-Out demo  \n " );
  printf( "------------------ \n" );
  printf( " * [u] -> Zoom in  \n" );
  printf( " * [d] -> Zoom out \n" );
  printf( " * [ESC] -> Close program \n \n" );

  /// Test image - Make sure it s divisible by 2^{n}
  string filename;
  filename = "/home/michael/pedestrian_detection/misc/custom_images/I00172_changed_equalized.png";
  filename = "/home/michael/pedestrian_datasets/caltech/data-INRIA/videos/set02/V000/00000002a.png";
  filename = "/home/michael/pedestrian_datasets/caltech/data-INRIA/videos/set01/V000/I00000.png";
  src = imread( filename );

  if( !src.data )
    { printf(" No data! -- Exiting the program \n");
      return -1; }

  tmp = src;
  dst = tmp;

  /// Create window
  namedWindow( window_name, CV_WINDOW_AUTOSIZE );
  imshow( window_name, dst );

  float ss = 1.00;
  float sr = 1.0718;

  double duration;
  duration = static_cast<double>(getTickCount());
  vector< pair<Mat,float> > scaleSpace;
  mk::getScaleSpace(src,scaleSpace,sr,ss);

  duration = static_cast<double>(getTickCount())-duration;
  duration /= getTickFrequency();
  printf( "** Scalespace created: time=%f, scales=%d \n",duration,scaleSpace.size());


  int curScale = 0;

  /// Loop
  while( true )
  {
    int c;
    c = waitKey(10);

    if( (char)c == 27 )
      { break; }
    if( (char)c == 'd' )
      {
    	curScale++;
        curScale = std::min((int)scaleSpace.size()-1,curScale);
    	printf( "** Scale: %d\n", curScale);

      }
    else if( (char)c == 'u' )
     {
       curScale--;
       curScale = std::max(0,curScale);
       printf( "** Scale: %d\n", curScale);
     }

    imshow( window_name, scaleSpace[curScale].first );
  }
  return 0;
}
