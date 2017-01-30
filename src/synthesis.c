
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdarg.h>
#include <string.h>

#include "quartic.h"

static double maxarg1,maxarg2;
#define MAX(a,b) (maxarg1=(a),maxarg2=(b),(maxarg1) > (maxarg2) ?\
	(maxarg1) : (maxarg2))

#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

#define MAX_BLOCK_LENGTH  1000
#define FFTDIV2           512
#define DELTA             80
#define OUTINC            4000
#define EXCINC            1000

short   count,debuglevel;

FILE *fp1,*fp2,*fp3,*fp4,*fp5,*fp6;

char timefile[100],pitchfile[100];
/****************************************************************************/

void display_matrix(double **a, int a_rows,int a_cols,char c[])
{
	int i,j;
		printf("%s\n",c);
		for (i=1;i<=a_rows;i++)
		{
			for (j=1;j<=a_cols;j++)	printf("%f\t",a[i][j]);
			printf("\n");
		}
		printf("\n");
}

/****************************************************************************/
void *MAT_malloc(size_t size)
{
  void *ptr;
if (size!=0)
{
  if ((ptr=malloc(size))==NULL)
  {
	 printf("\nOut Of Memory in MAT_malloc\n");

	 exit(0);
  }
}
else ptr=NULL;
  return(ptr);
}

/****************************************************************************/
void *MAT_realloc(void* ptr,size_t size)
{

  if ((ptr=realloc(ptr,size))==NULL)
  {
	 printf("\nOut Of Memory\n");
     
	 exit(0);
  }
  return(ptr);
}
/****************************************************************************/

double pythag(double a, double b)
{
	double absa,absb;
	absa=fabs(a);
	absb=fabs(b);
	if (absa>absb) return absa*sqrt(1.0+((absb/absa)*(absb/absa)));
	else return (absb=0.0 ? 0.0 : absb*sqrt(1.0+((absa/absb)*(absa/absb))));
}
/****************************************************************************/

void error(char error_text[])
{
	void exit();

	printf("Error...\n");
	printf("%s\n",error_text);
	printf("...press a key...\n");
	 
	exit(1);
}
/****************************************************************************/

double *vector(nl,nh)
int nl,nh;
{
	double *v;

	v=(double *)malloc((unsigned) (nh-nl+1)*sizeof(double));
	if (!v) error("allocation failure in vector()");
	return v-nl;
}


/****************************************************************************/

void free_vector(double *v,int nl,int nh)
{
	free((char*) (v+nl));
}
/****************************************************************************/
void free_matrix(double **m,int nrh,int nch)
{
	int nrl=0,ncl=0,i;
	
	for(i=nrh;i>=nrl;i--) free((char*) (m[i]+ncl));
	free((char*) (m+nrl));
}
/****************************************************************************/
double **matrix(nrh,nch)
int nrh,nch;
{
	int i,nrl=0,ncl=0;
	double **m;
	
	nrh++;
	nch++;
	
	m=(double **) malloc((unsigned) (nrh-nrl+1)*sizeof(double*));
	if (!m)
		error("allocation failure 1 in matrix()");

	m -= nrl;

	for(i=nrl;i<=nrh;i++)
	{
		m[i]=(double *) malloc((unsigned) (nch-ncl+1)*sizeof(double));
		if (!m[i]) 
			error("allocation failure 2 in matrix()");
		m[i] -= ncl;
	}
	return m;
}

/****************************************************************************/

void matrix_multiply
				(double** a,int a_rows,int a_cols,double **b,int b_rows,int b_cols,double** c)
{
	int i,j,z,c_rows,c_cols;
	
	if (a_cols!=b_rows)
		error("matrices must be compatible");		
	
	c_rows=a_rows;
	c_cols=b_cols;
		
	for(i=1;i<=c_rows;i++)
	{
		for(j=1;j<=c_cols;j++)
		{
			c[i][j]=0;
			for (z=1;z<=a_cols;z++)
			{
				c[i][j]+=a[i][z]*b[z][j];		
			}
		}
	}
} 
/****************************************************************************/
void transpose(double **a,int a_rows,int a_cols)
{
	double **c;
	int i,j;
	
	c=matrix(a_cols,a_rows);
	
	for (i=1;i<=a_rows;i++)
	{
		for(j=1;j<=a_cols;j++)
		{
			c[j][i]=a[i][j];			
		}
	}
	
	for (i=1;i<=a_rows;i++)
	{
		for(j=1;j<=a_cols;j++)
		{
			a[i][j]=c[i][j];			
		}
	}	

}
/****************************************************************************/

void svdcmp(double **a,int m,int n,double *w,double **v)

{
	int flag,i,its,j,jj,k,l,nm;
	double c,f,h,s,x,y,z;
	double anorm=0.0,g=0.0,scale=0.0;
	double *rv1,*vector();

	if (m < n) error("SVDCMP: You must augment A with extra zero rows");
	rv1=vector(1,n);
	for (i=1;i<=n;i++) {
		l=i+1;
		rv1[i]=scale*g;
		g=s=scale=0.0;
		if (i <= m) {
			for (k=i;k<=m;k++) scale += fabs(a[k][i]);
			if (scale) {
				for (k=i;k<=m;k++) {
					a[k][i] /= scale;
					s += a[k][i]*a[k][i];
				}
				f=a[i][i];
				g = -SIGN(sqrt(s),f);
				h=f*g-s;
				a[i][i]=f-g;
				if (i != n) {
					for (j=l;j<=n;j++) {
						for (s=0.0,k=i;k<=m;k++) s += a[k][i]*a[k][j];
						f=s/h;
						for (k=i;k<=m;k++) a[k][j] += f*a[k][i];
					}
				}
				for (k=i;k<=m;k++) a[k][i] *= scale;
			}
		}
		w[i]=scale*g;
		g=s=scale=0.0;
		if (i <= m && i != n) {
			for (k=l;k<=n;k++) scale += fabs(a[i][k]);
			if (scale) {
				for (k=l;k<=n;k++) {
					a[i][k] /= scale;
					s += a[i][k]*a[i][k];
			}
				f=a[i][l];
		g = -SIGN(sqrt(s),f);
				h=f*g-s;
				a[i][l]=f-g;
				for (k=l;k<=n;k++) rv1[k]=a[i][k]/h;
				if (i != m) {
					for (j=l;j<=m;j++) {
						for (s=0.0,k=l;k<=n;k++) s += a[j][k]*a[i][k];
						for (k=l;k<=n;k++) a[j][k] += s*rv1[k];
					}
				}
				for (k=l;k<=n;k++) a[i][k] *= scale;
			}
		}
		anorm=MAX(anorm,(fabs(w[i])+fabs(rv1[i])));
	}
	for (i=n;i>=1;i--) {
		if (i < n) {
			if (g) {
				for (j=l;j<=n;j++)
					v[j][i]=(a[i][j]/a[i][l])/g;
				for (j=l;j<=n;j++) {
					for (s=0.0,k=l;k<=n;k++) s += a[i][k]*v[k][j];
					for (k=l;k<=n;k++) v[k][j] += s*v[k][i];
				}
			}
			for (j=l;j<=n;j++) v[i][j]=v[j][i]=0.0;
		}
		v[i][i]=1.0;
		g=rv1[i];
		l=i;
	}
	for (i=n;i>=1;i--) {
		l=i+1;
		g=w[i];
		if (i < n)
			for (j=l;j<=n;j++) a[i][j]=0.0;
		if (g) {
			g=1.0/g;
			if (i != n) {
				for (j=l;j<=n;j++) {
					for (s=0.0,k=l;k<=m;k++) s += a[k][i]*a[k][j];
					f=(s/a[i][i])*g;
		for (k=i;k<=m;k++) a[k][j] += f*a[k][i];
				}
			}
			for (j=i;j<=m;j++) a[j][i] *= g;
	} else {
			for (j=i;j<=m;j++) a[j][i]=0.0;
		}
		++a[i][i];
	}
	for (k=n;k>=1;k--) {
		for (its=1;its<=30;its++) {
			flag=1;
			for (l=k;l>=1;l--) {
				nm=l-1;
				if ((double)(fabs(rv1[l])+anorm) == anorm) {
					flag=0;
					break;
				}
				if ((double)(fabs(w[nm])+anorm) == anorm) break;
			}
			if (flag) {
				c=0.0;
				s=1.0;
				for (i=l;i<=k;i++) {
					f=s*rv1[i];
					rv1[i]=c*rv1[i];
					if ((double)(fabs(f)+anorm) == anorm) break;
					g=w[i];
					h=pythag(f,g);
					w[i]=h;
					h=1.0/h;
					c=g*h;
					s=(-f*h);
					for (j=1;j<=m;j++) {
						y=a[j][nm];
						z=a[j][i];
						a[j][nm]=y*c+z*s;
						a[j][i]=z*c-y*s;
					}
				}
			}
			z=w[k];
			if (l == k) {
				if (z < 0.0) {
					w[k] = -z;
					for (j=1;j<=n;j++) v[j][k]=(-v[j][k]);
				}
				break;
			}
			if (its == 30) error("No convergence in 30 SVDCMP iterations");
		x=w[l];
			nm=k-1;
			y=w[nm];
			g=rv1[nm];
			h=rv1[k];
			f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
			g=pythag(f,1.0);
			f=((x-z)*(x+z)+h*((y/(f+SIGN(g,f)))-h))/x;
			c=s=1.0;
			for (j=l;j<=nm;j++) {
				i=j+1;
				g=rv1[i];
				y=w[i];
				h=s*g;
				g=c*g;
				z=pythag(f,h);
				
				rv1[j]=z;
				c=f/z;
				s=h/z;
				f=x*c+g*s;
				g=g*c-x*s;
				h=y*s;
				y=y*c;
				for (jj=1;jj<=n;jj++) {
					x=v[jj][j];
					z=v[jj][i];
					v[jj][j]=x*c+z*s;
					v[jj][i]=z*c-x*s;
				}
				z=pythag(f,h);
				w[j]=z;
				if (z) {
					z=1.0/z;
					c=f*z;
					s=h*z;
				}
				f=(c*g)+(s*y);
				x=(c*y)-(s*g);
				for (jj=1;jj<=m;jj++) {
					y=a[jj][j];
					z=a[jj][i];
					a[jj][j]=y*c+z*s;
					a[jj][i]=z*c-y*s;
				}
			}
		rv1[l]=0.0;
			rv1[k]=f;
			w[k]=x;
		}
	}
	free_vector(rv1,1,n);
}


/****************************************************************************/
short myRound(double x)
{
  short M;
  
  if (x>=0) M=(short)(x+0.5);
	 else    M=(short)(x-0.5);
    
  return(M);
}
/****************************************************************************/

double intcubic(double start,double lin,double quad,double cub,int n)
{
  return(start+lin*n+quad*n*n+cub*n*n*n);
}
/****************************************************************************/
double intlinear(double start,double lin,short n)
{
  return(start+lin*n);
}
/****************************************************************************/
void InitOutput(t_sentence *sent)
{
  sent->samples=MAT_malloc(OUTINC*sizeof(short));
  sent->spbufsize=OUTINC;
  sent->nsamps=0;

  sent->excits=MAT_malloc(EXCINC*sizeof(int));

  sent->exbufsize=EXCINC;
  sent->nexcits=0;
  sent->excits[0]=0;
}
/****************************************************************************/
void StoreExcits (t_sentence *sent,int *excits,short n)
{
  short i;

  if (sent->nexcits+n>sent->exbufsize)
  {
    sent->exbufsize+=EXCINC;
	 sent->excits=MAT_realloc(sent->excits,sent->exbufsize*sizeof(int));
  }


  for(i=0;i<n;i++)
    sent->excits[sent->nexcits+i]=excits[i];

  sent->nexcits+=n;
}
/****************************************************************************/
void StoreSamples(t_sentence *sent,double *speech,short n)
{
  short i;
  
  if (sent->nsamps+n>sent->spbufsize)
  {
    sent->spbufsize+=OUTINC;
    sent->samples=MAT_realloc(sent->samples,sent->spbufsize*sizeof(short));
  }

  for(i=0;i<n;i++)
	 sent->samples[sent->nsamps+i]=myRound(speech[i]);
    
  sent->nsamps+=n;
}
/****************************************************************************/
void FreeParams(t_params *params)
{
  int n;
  
  for (n=0;n<params->numframes;n++)
	 free(params->frame[n].peaks);

  free(params->frame);

}
/****************************************************************************/
void FreeSentence(t_sentence *sent)
{
  free(sent->samples);
  free(sent->excits );
}

/****************************************************************************/

void DumpToFile(t_sentence *sent)
{
  fwrite(sent->samples,sizeof(short),sent->nsamps,fp2);
  fclose (fp2);
  fwrite(sent->excits,sizeof(int),sent->nexcits,fp4);
  fclose(fp4);
}
/****************************************************************************/
void GetFrameLength(t_frame *Lframe,t_frame *Rframe,int *T,double t_scale)
{
  *T=(int)(t_scale*(Rframe->time-Lframe->time)+0.5);
}
/****************************************************************************/

double Average_Peak_Distance(double *position,short n)
{
  short i,sum=0;
  double F;
  
  for (i=1;i<n;i++)
    sum+=position[i]-position[i-1];
  
  F=sum/(double)(n-1);
  
  return(F);
}

/****************************************************************************/

void Delshort(short *position,short value,short num)

{
  short i;

  for (i=0;i<value;i++)     position[i]   = position[i];

  for (i=value+1;i<num;i++) position[i-1] = position[i];
}

/****************************************************************************/

void Match_Peaks(t_frame *Lframe,t_frame *Rframe,short *Lmatch,short *Rmatch)

{
  short Rcand=-1,
        Lcand=-1,
        Rcandpos=-1,
        Lcandpos=-1,

      N,M,n,m,n1,
      gap,smallest_gap,

      tentative,definitive,Lfreqi[MAXPEAKS],Rfreqi[MAXPEAKS];

  N=Lframe->numpeaks;
  M=Rframe->numpeaks;
  count=0;

  for (n=0;n<N;n++) Lfreqi[n]=(short)(Lframe->peaks[n].freq+0.5);
  for (m=0;m<M;m++) Rfreqi[m]=(short)(Rframe->peaks[m].freq+0.5);

  do
  {
	 n=0;

    smallest_gap=1000;

	 tentative=0;
    definitive=0;

    for (m=0;m<M;m++) /* SCAN K+1 FOR A CANDIDATE MATCH TO FIRST IN K */
	 {
      gap=abs(Lfreqi[n]-Rfreqi[m]);
      if ((gap<DELTA) && (gap<smallest_gap))
      {
        smallest_gap=gap;
		  Rcand =Rfreqi[m];
        Rcandpos  =m;
        tentative   =1;
      }
    }

    if (tentative==1) /* SCAN K FOR A DEFINITIVE MATCH TO K+1 */
    {
      for(n1=0;n1<N;n1++)
      {
        gap=abs(Lfreqi[n1]-Rfreqi[Rcandpos]);
		  if ((gap<DELTA) && (gap<smallest_gap+1))
		  {
          if ((gap==smallest_gap) && (n1>0))
          {

          }
          else
			 {
            smallest_gap=gap;
            Lcand =Lfreqi[n1];
            Lcandpos  =n1;
            definitive  =1;
          }
        }
		}
      if (Lcandpos !=0) /* IF THERE IS A BETTER MATCH TO K+1 IN K */
      {
        if ( (abs (Lfreqi[n]-Rfreqi[Rcandpos-1]) >=DELTA) && (Rcandpos>0) )
			 /* IF LOWER PEAK IN K+1 (IF ANY!!) IS MORE THAN DELTA AWAY....*/
        {
          tentative=0;  /* THIS SHOULD MAKE THE TRACK DIE..?*/
          definitive=0; /* WILL STOP IT DOING STUFF BELOW*/
        }
        if ((abs(Lfreqi[n]-Rfreqi[Rcandpos-1])<DELTA) && (Rcandpos>0) )
			 /* IF LOWER PEAK IN K+1 (IF ANY!!) IS LESS THAN DELTA AWAY....*/
		  {
			 Lcand =Lfreqi[n];
          Rcand =Rfreqi[Rcandpos-1];
          Lcandpos  =n;
          Rcandpos  =Rcandpos-1;
          definitive  =1;
        }

        if ((fabs(Lfreqi[n]-Rfreqi[Rcandpos-1])>=DELTA) && (Rcandpos==0) )
        {
          tentative=0;  /* THIS SHOULD MAKE THE TRACK DIE..?*/
          definitive=0; /* WILL STOP IT DOING STUFF BELOW*/
        }
      }
	 }
    if (definitive==1)    
    {
      Lmatch[count]=Rcand;
      Rmatch[count]=Lcand;
      count+=1;
      Delshort(&Lfreqi[0],Lcandpos,N);
      Delshort(&Rfreqi[0],Rcandpos,M);
		N=N-1;
      M=M-1;
    }

	 if (tentative==0)
    {
      Lmatch[count]=-1;
      Rmatch[count]=Lfreqi[0];
      count+=1;
		Delshort(&Lfreqi[0],0,N);
      N=N-1;
    }
  }
  while (N>0);

  for (m=0;m<M;m++)
  {
	 Lmatch[count]=Rfreqi[m];
    Rmatch[count]=-1;
    count+=1;
  }
  if ((Lframe->numpeaks==0)&&(Rframe->numpeaks==0)) count=0;
}

/****************************************************************************/

void Re_Order(t_frame *Lframe ,t_frame *Rframe,short *Lmatch,short *Rmatch)

/* This function converts the 'match' arrays from being a set of frequencies*/
/* to a set of indices......................................................*/

{
  short i,Lcount,Rcount,Lpeak,Rpeak;
  for (i=0;i<count;i++)
  {
    Lpeak=Rmatch[i];
    Rpeak=Lmatch[i];
    
    for (Lcount=0;Lcount<Lframe->numpeaks;Lcount++)
    {
		if ((short)(Lframe->peaks[Lcount].freq+0.5)==Lpeak)
      {
		  Rmatch[i]=Lcount;
        Lcount=Lframe->numpeaks;
      }
    }

	 for (Rcount=0;Rcount<Rframe->numpeaks;Rcount++)
    {
      if ((short)(Rframe->peaks[Rcount].freq+0.5)==Rpeak)
      {
        Lmatch[i]=Rcount;
        Rcount=Rframe->numpeaks;
		}
	 }
  }
}

/****************************************************************************/

void Generate_Speech(t_frame *Lframe,t_frame *Rframe,short *Lmatch,short *Rmatch,int T,int Z,double **inva,double *speech,short flag1,double gain)

{
  short   i;
  int 		Y,n,M;
  double  mag,phi,theta;
  double  *a1,*a2,*a3,*a4,*b1,*b2,**x,**y,*ampslope,m;
  t_pair w[MAXPEAKS],amp[MAXPEAKS],sysphase[MAXPEAKS],excphase[MAXPEAKS];

	double a,b,d,e,g,h,k,s;

  a1     =MAT_malloc(count*sizeof(double));
  a2     =MAT_malloc(count*sizeof(double));
  a3     =MAT_malloc(count*sizeof(double));
  a4     =MAT_malloc(count*sizeof(double));
  b1     =MAT_malloc(count*sizeof(double));
  b2     =MAT_malloc(count*sizeof(double));

	x=matrix(3,1);
  y=matrix(3,1);

  ampslope  =MAT_malloc(count*sizeof(double));

  for (i=0;i<count;i++)
  {
	 if (Rmatch[i]==-1)  /* TRACK IS TO BE BORN */
	 {
		w[i].R        =(2*PI/SAMPLERATE)*Rframe->peaks[Lmatch[i]].freq;
		w[i].L        =w[i].R;
		amp[i].R      =exp(Rframe->peaks[Lmatch[i]].amp);
		amp[i].L      =0.0;
		sysphase[i].R=Rframe->peaks[Lmatch[i]].sysphase;
		sysphase[i].L=sysphase[i].R;
			if(Rframe->type==0)
		  excphase[i].R=(double)rand()*(2*PI)/(double)RAND_MAX;
		else
			excphase[i].R	=0.0;
		excphase[i].L=intlinear(0.0,-w[i].R,Z);

	 }

	 if ((Lmatch[i]==-1)&&(Lframe->numpeaks>0))  /* TRACK IS TO BE KILLED */
	 {
		w[i].L        =(2*PI/SAMPLERATE)*Lframe->peaks[Rmatch[i]].freq;
		w[i].R        =w[i].L;
		amp[i].L      =exp(Lframe->peaks[Rmatch[i]].amp);
		amp[i].R      =0.0;
		sysphase[i].L=Lframe->peaks[Rmatch[i]].sysphase;
		sysphase[i].R=sysphase[i].L;
		excphase[i].L=Lframe->peaks[Rmatch[i]].excphase;
		excphase[i].R=intlinear(0.0,w[i].R,Z);
		if(Rframe->type==0)
		  excphase[i].R=(double)rand()*(2*PI)/(double)RAND_MAX;
	 }

	 if ((Lmatch[i]!=-1) && (Rmatch[i]!=-1)) /* TRACK IS MATCHED*/
	 {
		w[i].L        =(2*PI/SAMPLERATE)*Lframe->peaks[Rmatch[i]].freq;
		w[i].R        =(2*PI/SAMPLERATE)*Rframe->peaks[Lmatch[i]].freq;
		amp[i].L      =exp(Lframe->peaks[Rmatch[i]].amp);
		amp[i].R      =exp(Rframe->peaks[Lmatch[i]].amp);
		sysphase[i].L=Lframe->peaks[Rmatch[i]].sysphase;
		sysphase[i].R=Rframe->peaks[Lmatch[i]].sysphase;
		excphase[i].L=Lframe->peaks[Rmatch[i]].excphase;
		excphase[i].R=0.0;
		if(Rframe->type==0)
		  excphase[i].R=(double)rand()*(2*PI)/(double)RAND_MAX;
	 }


	 ampslope[i]=(amp[i].R-amp[i].L)/(double)T;

		a=inva[1][1];
		b=inva[1][2];
		d=inva[2][1];
		e=inva[2][2];
		g=inva[3][1];
		h=inva[3][2];

		if (Z>T) Y=Z; else Y=T;
   	
	//	k=a*a+3*Y*Y*d*d+(36/5.0)*Y*Y*Y*Y*g*g+3*Y*a*d+4*Y*Y*a*g+9*Y*Y*Y*d*g;

	//	s=2*a*b+6*Y*Y*d*e+(72/5.0)*Y*Y*Y*Y*g*h+3*Y*(a*e+b*d)+4*Y*Y*(a*h+b*g)+9*Y*Y*Y*(d*h+e*g);


	 m=(1/(2*PI))*(excphase[i].L+w[i].L*Z-excphase[i].R+(w[i].R-w[i].L)*Y*Y/(2.0*T));
	 M=myRound(m);

  //	 m=(1/(2*PI))*(excphase[i].L+w[i].L*Z-excphase[i].R-(w[i].R-w[i].L)*s/(2.0*k));
  //	 M=myRound(m);

		y[1][1]=excphase[i].R+2*PI*M-excphase[i].L-w[i].L*Z;
		y[2][1]=w[i].R-w[i].L;
		y[3][1]=0.0;

		matrix_multiply(inva,3,3,y,3,1,x);

		a1[i]  =x[1][1];
	 a2[i]  =x[2][1];
	 a3[i]  =x[3][1];

 /*	 a1[i]=(w[i].R-w[i].L)/(2*T);
	 a2[i]=0; // Quadratic interpolation
	 a3[i]=0;
   */
	 m=(-1.0/(2*PI))*(sysphase[i].R-sysphase[i].L);
	 M=myRound(m);

	 b1[i]  =(3.0/(T*T))    *(sysphase[i].R+2*PI*M-sysphase[i].L);
	 b2[i]  =(-2.0/(T*T*T)) *(sysphase[i].R+2*PI*M-sysphase[i].L);

  }

  for (n=0;n<T;n++)
  {
	 speech[n]=0.0;
	 for (i=0;i<count;i++)
	 {
		if (flag1)
		{
			mag   =intlinear(amp[i].L,ampslope[i],n);
		  phi   =sysphase[i].L+b1[i]*n*n+b2[i]*n*n*n;
		}
		else
		{
			mag=50;
			phi=0;
		}
		theta =excphase[i].L+w[i].L*n+a1[i]*n*n+a2[i]*n*n*n+a3[i]*n*n*n*n;
		if ((Lframe->type==0)&&(Rframe->type==0))
			theta+=((double)rand()*(0)/(double)RAND_MAX)-0;

		speech[n]+=mag*(cos(theta+phi));

	 }
  }

  for (i=0;i<count;i++)
  {

    if (Lmatch[i]!=-1) /* TRACK IS MATCHED OR BORN*/ 
		Rframe->peaks[Lmatch[i]].excphase=fmod((
      excphase[i].L+w[i].L*T+a1[i]*T*T+a2[i]*T*T*T+a3[i]*T*T*T*T),2*PI);
  }

  free(a1);
  free(a2);
  free(a3);
  free(a4);
  free(b1);
  free(b2);
  free(ampslope);
  free_matrix(x,3,1);
  free_matrix(y,3,1);

}
/****************************************************************************/
void Splint(double xa[],double ya[],double y2a[], short n, double x, double *y)

{
  short klo,khi,k;
  double h,b,a;
  
  klo=0;
  khi=n-1;
  while (khi-klo > 1)
  {
    k=(khi+klo) >> 1;
    if (xa[k] > x) khi=k;
    else klo=k;
  }

  h=xa[khi]-xa[klo];
  if (h==0.0) printf("Bad input to Splint");
  a=(double)(xa[khi]-x)/h;
  b=(double)(x-xa[klo])/h;
  
  *y=a*ya[klo]+b*ya[khi]+((a*a*a-a)*y2a[klo]+(b*b*b-b)*y2a[khi])*(h*h)/6.0;
}

/****************************************************************************/

void scale_pitch(t_frame *frame,double p_scale)
{
  short   i,num,oldnum;
  double  freq[MAXPEAKS],amp[MAXPEAKS],phase[MAXPEAKS],
  				mod_phase[MAXPEAKS],mod_freq[MAXPEAKS],mod_amp[MAXPEAKS],
  				y2[MAXPEAKS],y3[MAXPEAKS],bandwidth,lastgap;
    
  num=frame->numpeaks;
  oldnum=num;
  
  if (num>2) bandwidth=frame->peaks[num-1].freq;
  
  for(i=0;i<num;i++)
  {
	freq[i] =frame->peaks[i].freq;
		amp[i]  =frame->peaks[i].amp;
  	phase[i]=frame->peaks[i].sysphase; 
    y2[i]   =frame->peaks[i].y2;
    y3[i]   =frame->peaks[i].y3;
  }
  
  /* numpeaks = 2 + Actual Number */
  /* So go from 1 to N-1 */

  for (i=1;i<frame->numpeaks-1;i++)
  {
   	mod_freq[i-1]=p_scale*freq[i];      
	 Splint(freq,amp  ,y2,oldnum,mod_freq[i-1],&mod_amp[i-1]);
    Splint(freq,phase,y3,oldnum,mod_freq[i-1],&mod_phase[i-1]);
	if (mod_freq[i-1]>SAMPLERATE/2) num--;
  }
    
  num-=2; /* Remove the peaks at 0 and PI/2 */
    
  if ((num>2)&&(p_scale<1))	/* Add peaks to the end of the series */
  {
  	lastgap=mod_freq[num-1]-mod_freq[num-2];

	while(mod_freq[num-1]+lastgap<=bandwidth)
	{
  		mod_freq[num]=mod_freq[num-1]+lastgap;
  		Splint(freq,amp  ,y2,oldnum,mod_freq[num],&mod_amp[num]);
			Splint(freq,phase,y3,oldnum,mod_freq[num],&mod_phase[num]);
   	 	num++;
  	}
  }
  
  frame->numpeaks=num;   
  frame->peaks=realloc(frame->peaks,frame->numpeaks*sizeof(t_peak));
     
  for (i=0;i<frame->numpeaks;i++)
  {
	 frame->peaks[i].freq=mod_freq[i];
    frame->peaks[i].amp =mod_amp[i];
	 frame->peaks[i].sysphase=mod_phase[i];
  }
}
/****************************************************************************/
void GetData(t_params *params)
{
  short n;
  
  fread(&(params->numframes),sizeof(short),1,fp1);
  params->frame=MAT_malloc(params->numframes*sizeof(t_frame));
  for(n=0;n<params->numframes;n++)
  {
    fread(&(params->frame[n].numpeaks)  ,sizeof(short),1,fp1);
    fread(&(params->frame[n].time)      ,sizeof(int ),1,fp1);
    fread(&(params->frame[n].type)      ,sizeof(short),1,fp1);
	 fread(&(params->frame[n].period)    ,sizeof(short),1,fp1);
	 params->frame[n].peaks=MAT_malloc(params->frame[n].numpeaks*sizeof(t_peak));
    fread(params->frame[n].peaks,sizeof(t_peak),params->frame[n].numpeaks,fp1);
  }
}

/****************************************************************************/
void debug_pre(t_frame *frame,short i,double p_scale,double t_scale,int T, int Z)
{
    short n;
    
    fprintf(fp3,"Current        :%d\n",i);
    fprintf(fp3,"Pitch Scale    :%f\n",p_scale);
	 fprintf(fp3,"Time  Scale    :%f\n",t_scale);
    fprintf(fp3,"T              :%d\n",T);
    fprintf(fp3,"Z              :%d\n",Z);
    fprintf(fp3,"Type           :%d\n",frame->type);
    fprintf(fp3,"Original number:%d\n",frame->numpeaks);
    for (n=0;n<frame->numpeaks;n++)
	  fprintf(fp3,"Old freq:%f\tOld mag:%f\tOld phase:%f\n",frame->peaks[n].freq,frame->peaks[n].amp,frame->peaks[n].sysphase);
	 fprintf(fp3,"\n");
}

/****************************************************************************/
void debug_post(t_frame *frame)
{
	 short n;

    fprintf(fp3,"New number:%d\n",frame->numpeaks);
    for (n=0;n<frame->numpeaks;n++)
	  fprintf(fp3,"New freq:%f\tNew mag:%f\tNew phase:%f\n",frame->peaks[n].freq,frame->peaks[n].amp,frame->peaks[n].sysphase);
	 fprintf(fp3,"\n");
}
/****************************************************************************/
void GetExcits(t_sentence *sent,t_frame *Lframe,t_frame *Rframe,int *excits,short *nexcits,int T,int *Z,double p_scale)
{
 short i;
  static int lastexcit=0;

  int last,current,period;
  i=0;
  *nexcits=0;

  period=myRound((Rframe->time-Lframe->time)*(1.0/p_scale));

  last=lastexcit-sent->nsamps;

	if ((Rframe->type==0)&&(Lframe->type==0))
  {
		*nexcits=0;
	*Z=T;
  }

  if ((Lframe->type==0)&&(Rframe->type==1))
  {
	*nexcits=1;
	excits[0]=T;
	lastexcit=sent->nsamps+excits[*nexcits-1];
	*Z=T;
  }

  if ((Lframe->type==1)&&(Rframe->type==0))
  {
	*nexcits=0;
	*Z=T;
  }

  if ((Lframe->type==1)&&(Rframe->type==1))
  {
	 do
	 {
		current=last+period;
		excits[i]=current;
		i++;
		last=current;
	 }
	 while (current<=T);

		if(abs(current-T)<abs((current-period)-T)) *Z=current; else *Z=(current-period);

	  if (*Z==0) *Z=current;
	  if (Z<0) printf("\nZ neg\n");
	 *nexcits=i;

	 lastexcit=sent->nsamps+*Z;
  }
}
/****************************************************************************/
void Lagrange(int T,int Z,double **inva)
{
	//printf("In Lagrange: T:%d,Z:%d\n",T,Z);
	double **s,**p,**b,**invp,**invb,**a,**tmp,**v,w[5];
	int Y;
	short i,j;

	a=matrix(3,3);
	s=matrix(3,3);
	p=matrix(3,3);
	b=matrix(3,3);
	invp=matrix(3,3);
	invb=matrix(3,3);
	v=matrix(3,3);
    tmp=matrix(3,3);

	if(Z>T) Y=Z; else Y=T;

    a[1][1]=(double)(Z*Z); a[1][2]=(double)(Z*Z*Z);  a[1][3]=(double)(Z*Z)*(double)(Z*Z);
    a[2][1]=2*T;           a[2][2]=3*T*T;	         a[2][3]=4*T*T*T;
	a[3][1]=0.0;           a[3][2]=0.000;	         a[3][3]=0.00000;

	p[1][1]=1.0;           p[1][2]=1.5*Y;            p[1][3]=2.0*Y*Y;
	p[2][1]=0.0;           p[2][2]=sqrt(0.75)*Y;     p[2][3]=2*sqrt(0.75)*Y*Y;
	p[3][1]=0.0;           p[3][2]=0.0;              p[3][3]=sqrt(0.2)*Y*Y;    
		
	svdcmp(p,3,3,w,v);
			
	for(i=1;i<=3;i++)
	{
		for(j=1;j<=3;j++)
		{
			if(i==j) s[i][j]=1.0/w[i];
						else s[i][j]=0.0;
		}
	}
			
	transpose(p,3,3);
	matrix_multiply(s,3,3,p,3,3,tmp);     /* tmp=s.p (p is now u')*/
	matrix_multiply(v,3,3,tmp,3,3,invp);  /* invp=v.tmp*/ /* invp=v.s.u' */
			
	matrix_multiply(a,3,3,invp,3,3,b);

	svdcmp(b,3,3,w,v);
			
	for(i=1;i<=3;i++)
	{
		for(j=1;j<=3;j++)
		{
			if(i==j) s[i][j]=w[i];
						else s[i][j]=0;
		}	
	}
				
	for(i=1;i<=3;i++)
	{
		for(j=1;j<=3;j++)
		{
			if((i==j)&&(w[i]>0.0)) s[i][j]=1.0/w[i];
									 else s[i][j]=0.0;
		}
	}

	transpose(b,3,3);
	matrix_multiply(s,3,3,b,3,3,tmp);
	matrix_multiply(v,3,3,tmp,3,3,invb); /* invb=v.s.b' */
 	matrix_multiply(invp,3,3,invb,3,3,inva);

	inva[1][3]=0;
	inva[2][3]=0;
	inva[3][3]=0;

	free_matrix(a,3,3);
	free_matrix(s,3,3);
	free_matrix(p,3,3);
	free_matrix(b,3,3);
	free_matrix(invp,3,3);
	free_matrix(invb,3,3);
	free_matrix(v,3,3);
	free_matrix(tmp,3,3);

}
/****************************************************************************/
void OpenFiles()
{

 if((fp1=fopen(BINFILE,"rb"))==NULL)
	 error("Could not open BINFILE");

  fp3=fopen(DEBUGFILE,"w");
  fp4=fopen(EXFILE,"wb");


}
/****************************************************************************/
int main(int argc,char *argv[])
{
  short   i,Lmatch[MAXPEAKS],Rmatch[MAXPEAKS],nexcits,flag1,k,j;
  int    T,Z;
  double  *speech,*t_scale,*p_scale,**inva,power1,power2,gain;
  int    *excits;
  t_frame *Lframe;
  t_frame *Rframe;
	char filename[50];

  speech=MAT_malloc(MAX_BLOCK_LENGTH*sizeof(double));
  excits=MAT_malloc(10*sizeof(int));
  t_scale=MAT_malloc(100*sizeof(double));
  p_scale=MAT_malloc(100*sizeof(double));

	inva=matrix(3,3);

  if (argc!=4)
  {
	 printf("\nUsage: synthesis <time file> <pitch file> <speech/excitation?>\n");

	 exit(0);
  }

  strcpy(timefile,"c:\\work\\speech~1\\");
  strcat(timefile,argv[1]);
  strcat(timefile,".dat");

  strcpy(pitchfile,"c:\\work\\speech~1\\");
  strcat(pitchfile,argv[2]);
  strcat(pitchfile,".dat");

  strcpy(filename,"Q");
  strcat(filename,argv[1]);
  strcat(filename,"_");
  strcat(filename,argv[2]);
  strcat(filename,".pcm");

  flag1 =atoi(argv[3]);

  OpenFiles();

	if((fp5=fopen(timefile,"rb"))!=NULL)
	{
		fread(t_scale,sizeof(double),100,fp5);
		fclose(fp5);
	}
	else
	{
		printf("\nWarning: Time_scale constant (=%f)\n",atof(argv[1]));
		for (i=0;i<100;i++)
			t_scale[i]=atof(argv[1]);
	}

	if((fp6=fopen(pitchfile,"rb"))!=NULL)
	{
		fread(p_scale,sizeof(double),100,fp6);
		fclose(fp6);
	}
	else
	{
		printf("\nWarning: Pitch_scale constant (=%f)\n",atof(argv[2]));
		for (i=0;i<100;i++)
			p_scale[i]=atof(argv[2]);
	}


  InitOutput(&sent);
  GetData(&params);

/*  strcpy(filename,"quartic");
  strcat(filename,argv[1]);
  strcat(filename,"_");
  strcat(filename,argv[2]);
  if (flag1)
	strcat(filename,"_s.pcm");
  else
	strcat(filename,"_x.pcm");
									*/
  if ((fp2=fopen(filename,"wb"))==NULL) {printf("Cannot open output file %s",filename);exit(0);}

  printf("\nRead in %d Frames\n",params.numframes);
  printf("Please wait..Processing\n\n");


  for (i=0;i<params.numframes-1;i++)
  {
	 printf("Current Frame: %d\r",i);
	 fflush(stdout);

	 Lframe=&(params.frame[i]);
	 Rframe=&(params.frame[i+1]);

	 k=myRound(100*Lframe->time/(double)params.frame[params.numframes-1].time);
	 if (k<0) k=0;
	 if (k>99)k=99;

	 GetFrameLength(Lframe,Rframe,&T,t_scale[k]);

	 GetExcits(&sent,Lframe,Rframe,excits,&nexcits,T,&Z,p_scale[k]);

	 debug_pre(Rframe,i,p_scale[k],t_scale[k],T,Z);

	 power1=0;
	 for (j=0;j<Rframe->numpeaks;j++)
		power1+=exp(Rframe->peaks[j].amp)*exp(Rframe->peaks[j].amp);

	 if (Rframe->numpeaks>0)
	 {
		power1=power1/Rframe->numpeaks;
		power1=sqrt(power1);
	 }
	if ((Rframe->type==1)||(Lframe->type==1)||((i<params.numframes-2)&&(params.frame[i+2].type==1)))
	 scale_pitch(Rframe,p_scale[k]);
	 else
	 scale_pitch(Rframe,1.0);

	 power2=0;
	 for (j=0;j<Rframe->numpeaks;j++)
		power2+=exp(Rframe->peaks[j].amp)*exp(Rframe->peaks[j].amp);

	 if (Rframe->numpeaks>0)
	 {
		power2=power2/Rframe->numpeaks;
		power2=sqrt(power2);
	 }

	 if ((power1==0)||(power2==0))
	 gain=1;
	 else
	 gain=power1/power2;

    /**************************/
	 gain=1;

	 debug_post(Rframe);

	 Match_Peaks(Lframe,Rframe,Lmatch,Rmatch);

	 Re_Order   (Lframe,Rframe,Lmatch,Rmatch);

	 Lagrange(T,Z,inva);

	 Generate_Speech(Lframe,Rframe,Lmatch,Rmatch,T,Z,inva,speech,flag1,gain);

	 StoreSamples(&sent,speech,T);
	 StoreExcits (&sent,excits,nexcits);
  }
  printf("\nProcess Complete, Dumping...");
  fflush(stdout);
  DumpToFile(&sent);
  printf("ok\n");

  fclose(fp1);
  fclose(fp3);
  FreeParams(&params);
  FreeSentence(&sent);
  free (speech);
  free (excits);
  free_matrix(inva,3,3);

  }


