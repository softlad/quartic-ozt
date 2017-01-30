
#define   FFTPOINT         1024

typedef struct{
	
	double real;
	double imag;
	
}t_complex;

void SetupFFTInput(double* input,t_complex *fftdata,short N);
void FFT(t_complex *fftdata);
