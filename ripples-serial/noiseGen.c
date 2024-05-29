#include "omislib.h"

/*
*  Flip an array.
*/
void flip(_Dcomplex* a, int size) {
	_Dcomplex tmp;
	for (int i = 0; i < size/2; i++) {
		tmp = a[i];
		a[i] = a[size - 1 - i];
		a[size - 1 - i] = tmp;
	}
}

/*
* Calculate the inverse fourier transform.
* 
* @param *f: the signal array in frequency domain.
* @param *t: the result of the inverse fourier transform.
*/
void ifft(double *t, _Dcomplex *f, int size) {
	for (int n = 0; n < size; n++) {
		t[n] = 0;
		for (int k = 0; k < size; k++) {
			t[n] += creal(_Cmulcc(f[k], cexp(_Cbuild(0, (2*PI*n*k)/size)) ));
		}
		t[n] /= size;
	}
}

void noiseGen(double* noise, int size, double T, double dt) {
	int ff = 100; // [=] Hz, filter frequency
	double Tms = T * 1000;
	int s = ceil(Tms / dt);

	double dt_ins = dt / 1000;
	double df = 1 / (T+dt_ins); // frequency resolution
	double* fidx;
	double* faxis;
	double* filterf;

	double* Rr;
	int s1 = ceil(((float)s / 2));
	fidx = (double *)malloc(s1 * sizeof(double));
	faxis = (double*)malloc(s1 * sizeof(double));
	filterf = (double*)malloc(s1 * sizeof(double));
	_Dcomplex* fourier;
	fourier = (double*)malloc((s1 * 2 + 1) * sizeof(_Dcomplex));
	fourier[0] = _Cbuild(0, 0);

	_Dcomplex* fourierA;
	_Dcomplex* fourierB;
	fourierA = &(fourier[1]);
	fourierB = &(fourier[1+s1]);
	fourierA = (double*)malloc(s1 * sizeof(_Dcomplex));
	fourierB = (double*)malloc(s1 * sizeof(_Dcomplex));


	_Dcomplex phase_dist;

	Rr = malloc(s1 * sizeof(double));
	for (int i = 0; i < s1; i++) {
		fidx[i] = (double)i + 1;
		faxis[i] = (fidx[i] - 1) * df;

		Rr[i] = RandNormal();
		phase_dist = cexp(_Cbuild(0, PI * Rr[i]));
		filterf[i] = sqrt(1 / ( pow(2 * PI * ff,2) + pow(2 * PI * faxis[i],2) ));

		fourierA[i] = _Cmulcc(phase_dist, _Cbuild(filterf[i],0));
		fourierB[i] = conj(fourierA[i]);
	}
	flip(fourierB, s1);

	double* signal;
	signal = malloc((s1 * 2 + 1) * sizeof(double));
	ifft(signal, fourier, (s1 * 2 + 1));

	//print to check
	printf("%lf ", signal[1]);

	free(fidx);
	free(faxis);
	free(filterf);
	free(fourierA);
	free(fourierB);
	exit(0);
}