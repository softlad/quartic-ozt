package com.matpollard.speech;

public class Peak {	
	double freq;     /*Position of a peak (0-FFTPOINT/2-1) */
	double sysphase; /*Unwrapped System Phase of a peak    */
	double excphase; /*Excitation Phase of a peak          */
	double amp;      /*Amplitude of a peak                 */
	double y2;       /*First Derivitive for amplitude      */
	double y3;       /*First Derivitive for phase          */
	double deriv;    /*Second Derivitive of phase          */        
}