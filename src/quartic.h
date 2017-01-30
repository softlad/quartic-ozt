
#define TEXTFILE	"fcp/fcp.txt"
#define BINFILE	"fcp/fcp.bin"
#define EXFILE    "fcp/fcp.exc"
#define UPDATEFILE "fcp/fcp.upd"
#define DEBUGFILE "fcp/fcp.dbg"


#define SAMPLERATE        8000
#define DELTA             80
#define MAXPEAKS          1000
#define CUTOFF				  PI
#define BROWN				  0

typedef struct{

	short A;
	short E;
	short S;

}t_flag;

typedef struct{
	
	double *Acoefs;		/* Array of Amplitude Coefficients    */
	double *Ecoefs;		/* Array of Excit Phase Coefficients  */
	double *Scoefs;		/* Array of System Phase Coefficients */
	
	short		Aorder;		/* Orders for the above */
	short		Eorder;		/* NOTE the number of coeffs = 1 + order */
	short		Sorder;

}t_track;

typedef struct{
    
	int  beg;
	int  end;
	short avperiod;        
	short begindx;
	short endindx;

}t_section;

typedef struct{
	
	double freq;     /*Position of a peak (0-FFTPOINT/2-1) */
	double sysphase; /*Unwrapped System Phase of a peak    */
	double excphase; /*Excitation Phase of a peak          */
	double amp;      /*Amplitude of a peak                 */
	double y2;       /*First Derivitive for amplitude      */
	double y3;       /*First Derivitive for phase          */
	double deriv;    /*Second Derivitive of phase          */        

}t_peak;

typedef struct{
	
	t_peak *peaks;     /* An array of peaks              */
	short  numpeaks;   /* Number of peaks for a frame    */
	int   time;       /* Time at which analysis occured */
	short  type;       /* Voiced or Unvoiced             */
	short  period;		 /* Instantaneous Pitch Period     */
	 
}t_frame;
	
typedef struct{
	
	t_frame *frame;    /* An array of frames        */
	short   numframes; /* Number of frames analysed */
	short framebufsize;  /* Maximum frame buffer size */
	short nsections;
  short sectbufsize;
  t_section *section;

}t_params; 
	
t_params params;

typedef struct{
	
	short *samples;
	int  nsamps;
	int  spbufsize;
	int  *excits;
	int  nexcits;
	int  exbufsize;
	
}t_sentence;
	
t_sentence sent;
	


typedef struct{
	
	double L;
	double R;

}t_pair;
