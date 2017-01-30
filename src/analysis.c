
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "quartic.h"
#include "mathutils.h"
#include "memutils.h"
#include "dsp.h"

#define   FFTDIV2           512
#define   FRAMEINC          200
#define   SPEECHINC         8000   /* Enough for 1 second */
#define   EXCITSINC         200
#define   SECTINC           100
#define   TIMEGAP           0.005
#define   G                 3

FILE *fp1,*fp2,*fp3,*fp4,*fp5;
char *speech_file,*excits_file;

/****************************************************************************/
void InitParams (t_params *params)
{
  params->frame         = MAT_malloc(FRAMEINC*sizeof(t_frame));
  params->framebufsize  = FRAMEINC;
  params->numframes     = 0;
  params->nsections=0;
  params->sectbufsize=SECTINC;
  params->section   = MAT_malloc(SECTINC*sizeof(t_section));
}

/****************************************************************************/
void InitSentence (t_sentence *sent)
{
  sent->nsamps=0;
  sent->nexcits=0;
  sent->spbufsize=SPEECHINC;
  sent->exbufsize=EXCITSINC;
  sent->samples   = MAT_malloc(SPEECHINC*sizeof(short));
  sent->excits    = MAT_malloc(EXCITSINC*sizeof(int ));
}
/****************************************************************************/
void FreeParams(t_params *params)
{
  int n;
  
  for (n=0;n<params->numframes;n++)
    free(params->frame[n].peaks);

  free(params->frame);
  free(params->section);
}
/****************************************************************************/
void FreeSentence(t_sentence *sent)
{
  free(sent->samples);
  free(sent->excits );
}
/****************************************************************************/
void LoadSentence(t_sentence *sent)
{
  int num; 
  
  do
  {
    num=fread(&(sent->samples[sent->nsamps]),sizeof(short),SPEECHINC,fp2);
    if (num==SPEECHINC)
	 {
      sent->spbufsize+=SPEECHINC;
		sent->samples=MAT_realloc(sent->samples,sent->spbufsize*sizeof(short));
    }
    sent->nsamps+=num;
  }
  while (num==SPEECHINC);   
  
  do
  {
    num=fread(&(sent->excits[sent->nexcits]),sizeof(int),EXCITSINC,fp3);
    if (num==EXCITSINC)
    {
      sent->exbufsize+=EXCITSINC;
      sent->excits=MAT_realloc(sent->excits,sent->exbufsize*sizeof(int));
    }
    sent->nexcits+=num;
  }
  while (num==EXCITSINC);   

}

/****************************************************************************/

void StoreParams (t_frame *temp,t_params *params)
{
  short i;
  
  if (params->numframes + 1 >= params->framebufsize) 
  {
    params->framebufsize += FRAMEINC;
    params->frame = MAT_realloc(params->frame,params->framebufsize*sizeof(t_frame));
  } 

  params->frame[params->numframes].peaks=MAT_malloc(temp->numpeaks*sizeof(t_peak));
  
  for (i=0;i<temp->numpeaks;i++)
  {
	 params->frame[params->numframes].peaks[i].freq =temp->peaks[i].freq;
    params->frame[params->numframes].peaks[i].amp  =temp->peaks[i].amp;
    params->frame[params->numframes].peaks[i].sysphase=temp->peaks[i].sysphase;
    params->frame[params->numframes].peaks[i].excphase=temp->peaks[i].excphase;
    params->frame[params->numframes].peaks[i].y2   =temp->peaks[i].y2;
    params->frame[params->numframes].peaks[i].y3   =temp->peaks[i].y3;
  }
  
  params->frame[params->numframes].numpeaks=temp->numpeaks;
  params->frame[params->numframes].time    =temp->time;  
  params->frame[params->numframes].type    =temp->type;  
  params->frame[params->numframes].period  =temp->period;  

  params->numframes++;

}

/****************************************************************************/
void SetupWindow(double *input,double *output,short N)

/* This function sets up the non zero part of the fft analysis window */
{
  short   n;
  double *wind,sum;
  
  wind=(double*)MAT_malloc(N*sizeof(double));
  
  sum=0;

  for (n=0;n<N;n++)
  {
    wind[n]=0.5*(1-cos(2*PI*(n+1)/(1.0+N)));

    sum+=wind[n];
  }

  for (n=0;n<N;n++)
  {
	 wind[n]/=sum;
    output[n]=input[n]*wind[n];
  }
  free(wind);
}
/****************************************************************************/
short WhatSection(t_params *params,int time)
{
  short i,section;
  
  if (time<0)
  {
    printf("error in What Section\n");
    exit(0);
  }
  
  section=-1;
  
  for (i=0;i<params->nsections;i++)
	 if ((params->section[i].beg<=time) && (params->section[i].end>=time)) section=i;
    
  return(section);
}

/****************************************************************************/
void SplinePhases(t_frame *temp)

{
  short i,k,n;
  double p,qn,sig,un,*u,yp1,ypn;
  
  n=temp->numpeaks;
    
  u=MAT_malloc(n*sizeof(double));
    
  yp1=0.0;          /* These are definable */
  ypn=1e30;
  
  if (yp1 >0.99e30)
	 temp->peaks[0].y3=u[0]=0.0;
  else
  {
    temp->peaks[0].y3=-0.5;
    u[0]=(3.0/(temp->peaks[1].freq-temp->peaks[0].freq))*((temp->peaks[1].sysphase-temp->peaks[0].sysphase)/(temp->peaks[1].freq-temp->peaks[0].freq)-yp1);
  }
  
  for (i=1;i<n-1;i++)
  {
    sig=(double)(temp->peaks[i].freq-temp->peaks[i-1].freq)/(temp->peaks[i+1].freq-temp->peaks[i-1].freq);
    p=sig*temp->peaks[i-1].y3+2.0;
    temp->peaks[i].y3=(sig-1.0)/p;
    u[i]=(temp->peaks[i+1].sysphase-temp->peaks[i].sysphase)/(temp->peaks[i+1].freq-temp->peaks[i].freq)
    -(temp->peaks[i].sysphase-temp->peaks[i-1].sysphase)/(temp->peaks[i].freq-temp->peaks[i-1].freq);
    u[i]=(6.0*u[i]/(temp->peaks[i+1].freq-temp->peaks[i-1].freq)-sig*u[i-1])/p;
  }
  
  if (ypn > 0.99e30)
    qn=un=0.0;
  else
  {
    qn=0.5;
    un=(3.0/(temp->peaks[n-1].freq-temp->peaks[n-2].freq))* (ypn-(temp->peaks[n-1].sysphase-temp->peaks[n-2].sysphase)/(temp->peaks[n-1].freq-temp->peaks[n-2].freq));
  }
  
  temp->peaks[n-1].y3=(un-qn*u[n-2])/(qn*temp->peaks[n-2].y3+1.0);
  for (k=n-2;k>=0;k--)
    temp->peaks[k].y3=temp->peaks[k].y3*temp->peaks[k+1].y3+u[k];
    
  free(u);
}

/****************************************************************************/
void SplineAmps(t_frame *temp)

{
  short i,k,n;
  double p,qn,sig,un,*u,yp1,ypn;

  n=temp->numpeaks;
    
  u=MAT_malloc(n*sizeof(double));

    
  yp1=0.0;          /* These are definable */
  ypn=1e30;
  
  if (yp1 >0.99e30)
    temp->peaks[0].y2=u[0]=0.0;
  else
  {
    temp->peaks[0].y2=-0.5;
    u[0]=(3.0/(temp->peaks[1].freq-temp->peaks[0].freq))*((temp->peaks[1].amp-temp->peaks[0].amp)/(temp->peaks[1].freq-temp->peaks[0].freq)-yp1);
  }
  
  for (i=1;i<n-1;i++)
  {
	 sig=(double)(temp->peaks[i].freq-temp->peaks[i-1].freq)/(temp->peaks[i+1].freq-temp->peaks[i-1].freq);
    p=sig*temp->peaks[i-1].y2+2.0;
    temp->peaks[i].y2=(sig-1.0)/p;
    u[i]=(temp->peaks[i+1].amp-temp->peaks[i].amp)/(temp->peaks[i+1].freq-temp->peaks[i].freq)
    -(temp->peaks[i].amp-temp->peaks[i-1].amp)/(temp->peaks[i].freq-temp->peaks[i-1].freq);
    u[i]=(6.0*u[i]/(temp->peaks[i+1].freq-temp->peaks[i-1].freq)-sig*u[i-1])/p;
  }
  
  if (ypn > 0.99e30)
    qn=un=0.0;
  else
  {
    qn=0.5;
    un=(3.0/(temp->peaks[n-1].freq-temp->peaks[n-2].freq))* (ypn-(temp->peaks[n-1].amp-temp->peaks[n-2].amp)/(temp->peaks[n-1].freq-temp->peaks[n-2].freq));
  }
  
  temp->peaks[n-1].y2=(un-qn*u[n-2])/(qn*temp->peaks[n-2].y2+1.0);
  for (k=n-2;k>=0;k--)
    temp->peaks[k].y2=temp->peaks[k].y2*temp->peaks[k+1].y2+u[k];

  free(u);
}

/****************************************************************************/

void PickPeaks(t_params *params,t_frame *temp,t_complex *fftdata)

{
  short   section,i,F,beg,end,last_position,found=0,pos_of_biggest_so_far;

  double  mag1,mag2,mag3,amp[FFTDIV2],pha[FFTDIV2],biggest_so_far,
          phase_of_biggest_so_far;
  
  temp->numpeaks=1;
  
  section=WhatSection(params,temp->time);
  if (section!=-1) 
    F=myRound(FFTPOINT/(float)params->section[section].avperiod);
  else 
    F=1;
         
  atan4(0,0,1);
  
  for(i=0; i<FFTDIV2; i++)
  {
    if ((fftdata[i].real==0)&&(fftdata[i].imag==0))
		amp[i]=0;
    else
      amp[i]=log(sqrt(sqr(fftdata[i].real)+sqr(fftdata[i].imag)));
    
    if ((fftdata[i].imag!=0)&&(fftdata[i].real!=0)) 
      pha[i]=atan4(fftdata[i].imag,fftdata[i].real,0);        
    else pha[i]=0;
  }
      
  last_position=0;
  while((last_position+F+(F/2))<FFTDIV2)
  {
    beg=(short)( (last_position+F-F/2)+0.5);
    end=(short)( (last_position+F+F/2)+0.5);

    found=0;
    biggest_so_far=0;
    for(i=beg;i<=end;i++)
    {
		mag1=amp[i-1];
      mag2=amp[i];
      mag3=amp[i+1];
          
      if ((mag1<mag2) && (mag3<mag2) && (mag2>biggest_so_far))
      {         
        pos_of_biggest_so_far=i;
        biggest_so_far=mag2;
        phase_of_biggest_so_far=pha[i];
        found=1;
      }   
    }
    
    if (found==1)
    {
      temp->peaks[temp->numpeaks].freq  =SAMPLERATE*pos_of_biggest_so_far/(double)FFTPOINT;
      temp->peaks[temp->numpeaks].amp   =biggest_so_far;
      temp->peaks[temp->numpeaks].sysphase =phase_of_biggest_so_far;
      temp->peaks[temp->numpeaks].excphase =0;
		temp->numpeaks++;
      last_position=pos_of_biggest_so_far;
    }     
    if (found==0) last_position+=F;
  }
  
   temp->peaks[0].freq=0;
   
   if(temp->numpeaks==1)
   	temp->peaks[0].amp=amp[0];
   else
   	temp->peaks[0].amp=amp[0]+(temp->peaks[1].amp-amp[0])/2.0;
   
   temp->peaks[0].sysphase=temp->peaks[1].sysphase;
   temp->peaks[0].excphase=0;
   
   temp->peaks[temp->numpeaks].freq=SAMPLERATE/2.0;
   if(temp->numpeaks==1)
   	temp->peaks[temp->numpeaks].amp=amp[FFTDIV2-1];
	else
   	temp->peaks[temp->numpeaks].amp=temp->peaks[temp->numpeaks-1].amp+
   	(amp[FFTDIV2-1]-temp->peaks[temp->numpeaks-1].amp)/2.0;
   temp->peaks[temp->numpeaks].sysphase=temp->peaks[temp->numpeaks-1].sysphase;
   temp->peaks[temp->numpeaks].excphase=0;
   temp->numpeaks++;
   
}

/****************************************************************************/
void GetBlockLength(t_frame *temp,t_params *params,short *N)
{ 
  short section,avperiod,sum=0;
  
  section=WhatSection(params,temp->time);
  
  if (section!=-1)
  {
  	avperiod=params->section[section].avperiod;
  }
  else
  {
  	for (section=0;section<params->nsections;section++) 
  		sum+=params->section[section].avperiod;
  	avperiod=myRound((double)sum/params->nsections);
  }
  
  *N=myRound(2.5*avperiod);
  if ((*N)%2==0) (*N)++; 
}
/****************************************************************************/
void GetNextTime(t_frame *temp,t_sentence *sent,short updategap)
{
	short i=-1,type,prev_voiced=0;
	int gap1,gap2,Z1,Z2,t;
	static int lasttime=0;
	
		t=temp->time;
		/* Find Nearest Excit Greater than current time */
		do 
		{
			i++; 
			Z1=sent->excits[i];
		}
		while(Z1<=t);

		if(i==0) prev_voiced=0;
		else
			if(sent->excits[i-1]==t) prev_voiced=1;
		
		gap1=Z1-t;
		
		if(i==sent->nexcits-1) gap2=0;
		else
		{
			Z2=sent->excits[i+1];
			gap2=Z2-Z1;
		}
	
		if(gap1>2*updategap) type=1;
		if( (!prev_voiced) && (gap1>updategap) && (gap1<=2*updategap) && (gap2<=2*updategap) ) type=2;
		if( (!prev_voiced) && (gap1>updategap) && (gap1<=2*updategap) && (gap2>2*updategap)  ) type=1;
		if( (prev_voiced)  && (gap1>updategap) && (gap1<=2*updategap)) type=3;
		if( (prev_voiced)  && (gap1<=updategap)) type=3;
		if( (!prev_voiced) && (gap1<=updategap) && (gap2<=2*updategap) ) type=3;
		if( (!prev_voiced) && (gap1<=updategap) && (gap2>2*updategap) ) type=1;
	
		temp->type=0;
		temp->period=0;
		
		switch(type)
		{	
			case 1:	temp->time+=updategap;
							if ((gap1==updategap)&&(gap2<=2*updategap))
							{
								temp->type=1;
								temp->period=gap2;
							}
							break;
							
			case 2: temp->time=sent->excits[i]-(sent->excits[i]-temp->time)/2;
							break;
							 
			case 3: temp->time=sent->excits[i];
							if (gap2<=2*updategap)
							{
								temp->period=myRound((Z1-t)+(Z2-2*Z1+t)*0.5*(Z1-t)/( (Z2-0.5*(Z2-Z1))-(Z1-0.5*(Z1-t)) ));
								temp->type=1;
							}
							else
							{
								temp->period=gap1;
								temp->type=1;
							}
							break;

			default: printf("Error in assignment of next update point\n");
								exit(0);
		}	
		if ((lasttime<sent->nsamps-1)&&(temp->time>sent->nsamps-1)) 
			temp->time=sent->nsamps-1;
		
		lasttime=temp->time;
}
/****************************************************************************/
void GetUpdateGap(t_params *params,short *updategap)
{
	short sum=0,n=0,i;
	
	for(i=0;i<params->nsections;i++)
	{
		if(params->section[i].avperiod>0) 
		{
			sum+=params->section[i].avperiod;
			n++;
		}
	}
	*updategap=myRound(sum/(float)n);
}
/****************************************************************************/
void ReadData(t_frame *temp,t_sentence *sent,double *x,short N)
{ 
  short mid,r_ntoread,l_ntoread;
  int n;
     
  mid=(N-1)/2;
  r_ntoread=(N+1)/2;
  l_ntoread=(N-1)/2;
          
  for (n=temp->time;n<temp->time+r_ntoread;n++)
    x[n-temp->time+mid]=(double)sent->samples[n];
        
  if (temp->time>mid)
  {
    for (n=0;n<l_ntoread;n++)
      x[n]=(double)(sent->samples[temp->time-l_ntoread+n]); 
  }
  
  else
  {
   for (n=0;n<(mid-temp->time);n++)
      x[n]=0.0;
    for (n=0;n<temp->time;n++)
      x[n+(mid-temp->time)]=(double)sent->samples[n]; 
  } 
}
/****************************************************************************/
void OutputUpd(t_params *params)
{
 	int n;
 	
 	for(n=0;n<params->numframes;n++)
		fwrite(&(params->frame[n].time),sizeof(int),1,fp5);
}
/****************************************************************************/
void OutputBin(t_params *params)
{
  short n;
  
  fwrite(&(params->numframes),sizeof(short),1,fp4);
  
  for(n=0;n<params->numframes;n++)
  {
    fwrite(&(params->frame[n].numpeaks) ,sizeof(short),1,fp4);
    fwrite(&(params->frame[n].time)     ,sizeof(int ),1,fp4);
    fwrite(&(params->frame[n].type)     ,sizeof(short),1,fp4);
    fwrite(&(params->frame[n].period)   ,sizeof(short),1,fp4);
    
    fwrite(params->frame[n].peaks,sizeof(t_peak),params->frame[n].numpeaks,fp4);
  }
}
/****************************************************************************/
void OutputTxt(t_params *params)
{
  short n,i;
   
  fprintf(fp1,"\n\n");
  fprintf(fp1,"Test output for Analysis Module...\n");
  fprintf(fp1,"Number of Frames   :%d\n",params->numframes);
  fprintf(fp1,"Number of Sections :%d\n",params->nsections);
  fprintf(fp1,"\n");
  for(i=0;i<params->nsections;i++)
  {
    fprintf(fp1,"Section   %d\n",i);
    fprintf(fp1,"From      %d samples to %d samples\n",params->section[i].beg,params->section[i].end);
    fprintf(fp1,"Av Period %d\n",params->section[i].avperiod);
    fprintf(fp1,"\n");
  }

  for (n=0;n<params->numframes;n++)
  {
    fprintf(fp1,"\tFrame Number   :%d\n",n);
    fprintf(fp1,"\tTime (samples) :%d\n",params->frame[n].time);
    fprintf(fp1,"\tNumber of Peaks:%d\n",params->frame[n].numpeaks);
    fprintf(fp1,"\tVoiced?        :%d\n",params->frame[n].type);
    fprintf(fp1,"\tPeriod         :%d\n",params->frame[n].period);

	 for(i=0;i<params->frame[n].numpeaks;i++)
	 {
		fprintf(fp1,"\tFreq:%.0f\tAmp:%.4f\tSysPhase:%.4f\tExcPhase:%.4f\ty2:%.6f\ty3:%.6f\n", params->frame[n].peaks[i].freq,params->frame[n].peaks[i].amp,params->frame[n].peaks[i].sysphase,params->frame[n].peaks[i].excphase,params->frame[n].peaks[i].y2,params->frame[n].peaks[i].y3);
	 }
	 fprintf(fp1,"\n");

  }
}
/****************************************************************************/
void SectionSpeech(t_sentence *sent,t_params *params)
{
  int smallest_gap=sent->nsamps,gap,A,B,C,D,E;

  short i,sum,L1,L2,L3,L4,beg,end,in,count;

  //printf("Smallest Gap is %d\n", smallest_gap);

  params->nsections=0;

  /* Determine Smallest Gap Between Adjacent Excitation Points */
  for (i=0;i<sent->nexcits;i++)
  {
    //printf("%d -> ",sent->excits[i]);
    //printf("%d",sent->excits[i+1]);
    //printf("\n");

	 gap=sent->excits[i+1]-sent->excits[i];
   //printf("Gap is %d\n", gap);
	 if ((gap>=0)&&(gap<smallest_gap)&&(gap!=0)) smallest_gap=gap;
  }

  //printf("\nSmallest Gap is %d", smallest_gap);

  if(sent->nexcits==0)
  {
	 printf("Error: during SectionSpeech\n");
	 exit(0);
  }

  for (i=0;i<sent->nexcits;i++)
  {
	if (i<2) A=0; else A=sent->excits[i-2];
	if (i<1) B=0; else B=sent->excits[i-1];
	C=sent->excits[i];
	if (i>sent->nexcits-2) D=sent->nsamps;
			  else D=sent->excits[i+1];
	if (i>sent->nexcits-3) E=sent->nsamps;
			  else E=sent->excits[i+2];

	L1=((B-A)>G*smallest_gap);
	L2=((C-B)>G*smallest_gap);
	L3=((D-C)>G*smallest_gap);
	L4=((E-D)>G*smallest_gap);

	if (i==0) L1=L2=1;
	if (i==1) L1=1;
	if (i==sent->nexcits-2) L4=1;
	if (i==sent->nexcits-1) L3=L4=1;

		in=(((!L1)&&(!L2)) || ((!L3)&&(!L4)) || ((!L2)&&(!L3)) );
		beg=((L2)&&(!L3)&&(!L4));
		end=((!L1)&&(!L2)&&(L3));

		if (beg)
		{
			sum=0;
			count=0;
			params->section[params->nsections].beg=C;
			params->section[params->nsections].begindx=i;
		}

		if ((in)&&(!beg))
		{
			sum+=C-B;
			count++;
		}

		if (end)
		{
			params->section[params->nsections].end=C;
			params->section[params->nsections].endindx=i;
			params->section[params->nsections].avperiod=myRound((double)sum/count);
			params->nsections++;
		}
  }

	if (params->nsections==0)
  {
	printf("\nError..no sections found\n");
	exit(0);
  }
  params->section[params->nsections].end=sent->samples[sent->nsamps];
}

/****************************************************************************/
void OpenFiles()
{
  fp1=fopen(TEXTFILE,"w");
  if( (fp2=fopen(speech_file,"rb"))==NULL)
  {
	 printf("\nCannot open speech file\n");
	 //getch();
	 exit(0);
  }
  if((fp3=fopen(excits_file,"rb"))==NULL)
  {
    printf("\nCannot open .exc file\n");
    //getch();
	 exit(0);
  }
  fp4=fopen(BINFILE,"wb");
  fp5=fopen(UPDATEFILE,"wb");
}
/****************************************************************************/
int main(int argc,char *argv[])
{
	short			N,updategap;
	double		*x,*xfftwind;
	t_frame		temp;
	t_complex *fftdata;

	fftdata			=	(t_complex*)		MAT_malloc(FFTPOINT*sizeof(t_complex));
	x					=	(double*)			MAT_malloc(FFTPOINT*sizeof(double));
	xfftwind			=	(double*)			MAT_malloc(FFTPOINT*sizeof(double));
	temp.peaks		=	(t_peak*)			MAT_malloc(MAXPEAKS*sizeof(t_peak));

	if (argc!=3)
	{
		printf("\nUsage: fcp_analy <speechfile> <pitchmarks> \n");
      //getch();
		exit(0);
	}
	
	speech_file=argv[1];
	excits_file=argv[2];
	OpenFiles();

  //printf("Init Params\n");
	InitParams(&params);
  //printf("Init Sentence\n");
	InitSentence(&sent);
  //printf("Load Sentence\n");
	LoadSentence(&sent);

  printf("\nLoaded in %f Seconds of Speech..\n",sent.nsamps/(float)SAMPLERATE);
  printf(  "Loaded in %d Speech Samples.....\n",sent.nsamps);
  printf(  "Loaded in %d Excitation Points..\n",sent.nexcits);
  printf(  "Please Wait Processing...\n\n");



  //printf("Section Speech\n");
	SectionSpeech(&sent,&params);
  //printf("Get Update Gap\n");
	GetUpdateGap(&params,&updategap);


	temp.time=0;

	/* Beginning of main loop */
	do
	{
		printf("Completed %d%%\r",(int)(100*(float)temp.time/sent.nsamps));
		fflush(stdout);

		GetBlockLength(&temp,&params,&N);

		ReadData(&temp,&sent,x,N);

		SetupWindow(x,xfftwind,N);

		SetupFFTInput(xfftwind,fftdata,N);

		FFT(fftdata);

		PickPeaks(&params,&temp,fftdata);

		if(temp.numpeaks>1)
		{
			SplinePhases(&temp);
			SplineAmps(&temp);
		}

		StoreParams(&temp,&params);

		GetNextTime(&temp,&sent,updategap);
	}
	while (temp.time<sent.nsamps);

	printf("Completed 100%%\nDumping...");
	fflush(stdout);

	OutputTxt(&params);
	OutputBin(&params);
	OutputUpd(&params);
	printf("Ok\n");
	fflush(stdout);

	fclose(fp1);
	fclose(fp2);
	fclose(fp3);
	fclose(fp4);
	fclose(fp5);

	FreeParams(&params);
	FreeSentence(&sent);
	free(x);
	free(xfftwind);
	free(fftdata);
	free(temp.peaks);
}

