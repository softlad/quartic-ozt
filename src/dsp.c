#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mathutils.h"
#include "dsp.h"

/****************************************************************************/
void SetupFFTInput(double* input,t_complex *fftdata,short N)

/* This function sets up the entire fft analysis window */
/* Input is the windowed excitation signal              */
/* Output is the zero-padded fft input                  */

{
  short n,mid,next;
  
  mid=(N-1)/2;
   
  for (n=0;n<=(N-1)/2;n++)
  {
    fftdata[n].real=input[n+(N-1)/2];
   fftdata[n].imag=0.0;
  } 
  /* Place the right hand side of the window at the beginning */

  next=1+(N-1)/2;

  for (n=next;n<next+(FFTPOINT-N);n++)
  {
    fftdata[n].real=0.0;
    fftdata[n].imag=0.0;   
  }
  /* Place the left hand side of the window at the end        */
  
  next=next+(FFTPOINT-N);

  for (n=next;n<FFTPOINT;n++)
  {
     fftdata[n].real=input[(n-next)];
     fftdata[n].imag=0.0;   
  }
  /* Fill in between with zeros */

}

/****************************************************************************/
void FFT(t_complex *fftdata)

/* Input  is Real I/O and Imaginary I/O              */
/* Output is Real I/O and Imaginary I/O              */

{
  int N,N1,N2,NM,N2M,i,j,k,l;
  double e,c,s,a,xt,yt;

  N=FFTPOINT;      

  N2=N ;
  NM=N-1;
  k=2;
  e= PI/(float)N;
  do
  {
    N1=N2;
    N2=N2/2;
    N2M=N2-1 ;
    e=e*2.0;
    a=0.0;
    for(j=0;j< N2M+1; j++)
    {
      c=cos(a);
      s=sin(a);
      a=a+e;
      i=j;

      do
      {
        l=i+N2;
        xt=fftdata[i].real-fftdata[l].real;
        fftdata[i].real+=fftdata[l].real;

        yt=fftdata[i].imag-fftdata[l].imag;
        fftdata[i].imag+=fftdata[l].imag;
        
        fftdata[l].real=xt*c+yt*s;
        fftdata[l].imag=yt*c-xt*s;
        i=i+N1;
      }
      while(i < (NM+1));
    }
    k=k+k;  /*end of J loop */
  }
  while( k < (N+1));

  j=0;

  for(i=0; i<N-1; i++)
  {
    if(i<j)
	 {
      xt=fftdata[j].real;
      fftdata[j].real=fftdata[i].real;
      fftdata[i].real=xt;
      yt=fftdata[j].imag;
      fftdata[j].imag=fftdata[i].imag;
      fftdata[i].imag=yt;
    }
    k=N/2;
    while((k-1)<j) { j=j-k; k=k/2;}  j=j+k;
  }    /* end of i loop */

 }   /* end of procedure fft */