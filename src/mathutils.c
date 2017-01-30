#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mathutils.h"

/****************************************************************************/
short myRound(double x)
{
  short M;
  
  if (x>=0) M=(short)(x+0.5);
    else    M=(short)(x-0.5);
    
  return(M);
}
/***************************************************************************/
double sqr(double val)
{
  double ans;
  ans=val*val;

  return(ans);
}
/****************************************************************************/
double atan4(double im,double re,short reset)
{
  static short  i=0;
  static double Old=0.0;
  double        New,angle;

  if (reset) 
  {
    i=0;
    Old=0.0;
    return(0);
  }
  else
  {
    New=atan2(im,re);
	 if ((Old>PI/2)&&(Old<PI)&&(New>-PI)&&(New<-PI/2)) i++;
    if ((Old>-PI)&&(Old<-PI/2)&&(New>PI/2)&&(New<PI)) i--;

    angle=(2*i*PI)+New;
    Old=New;

    return(angle);
  }
}