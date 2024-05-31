#include "omislib.h"

/*
* Calculate the inverse fourier transform.
* 
* @param *f: the signal array in frequency domain.
* @param *t: the result of the inverse fourier transform.
*/
void idft(double *t, _Dcomplex *f, int size) {
	for (int n = 0; n < size; n++) {
		t[n] = 0;
		for (int k = 0; k < size; k++) {
			t[n] += creal(_Cmulcc(f[k], cexp(_Cbuild(0, (2*PI*n*k)/size)) )); //throwing away half of the info?
		}
		t[n] /= size;
	}
}

double calculateMean(int N, double *data)
{
	double sum = 0;
	for (int i = 0; i < N; i++) {    
		sum += data[i];
	}

	return sum / N;
}

double calculateStandardDeviation(int N, double* data)
{
	double sum = 0;
	for (int i = 0; i < N; i++) {
		sum += data[i];
	}

	double mean = sum / N;
	double values = 0;

	for (int i = 0; i < N; i++) {
		values += pow(data[i] - mean, 2);
	}
	double variance = values / N;

	return sqrt(variance);
}

void noiseGen(double* noise, int size, double T, double dt) { //size must be odd
	int ff = 100; // [=] Hz, filter frequency

	double df = 1 / (T+(dt / 1000)); // frequency resolution
	double* fidx;
	double* faxis;
	double* filterf;

	double* Rr;
	int s1 = (size-1) / 2 ;
	fidx = (double *)malloc(s1 * sizeof(double));
	faxis = (double*)malloc(s1 * sizeof(double));
	filterf = (double*)malloc(s1 * sizeof(double));
	_Dcomplex* fourier;
	fourier = (double*)malloc(size * sizeof(_Dcomplex));
	//fourier[s1] = _Cbuild(0, 0);
	fourier[0] = _Cbuild(0, 0);

	_Dcomplex phase_dist;
	
	_Dcomplex tmp;

	Rr = malloc(s1 * sizeof(double));
	for (int i = 0; i < s1; i++) {
		fidx[i] = (double)i + 1;
		faxis[i] = (fidx[i] - 1) * df;

		Rr[i] = RandNormal();
		phase_dist = cexp(_Cbuild(0, PI * Rr[i]));   
		filterf[i] = sqrt(1 / ( pow(2 * PI * ff,2) + pow(2 * PI * faxis[i],2) ));

		tmp = _Cmulcc(phase_dist, _Cbuild(filterf[i], 0));

		fourier[1 + i] = tmp;
		fourier[s1 + 1 + i] = conj(tmp);
	}

	idft(noise, fourier, size);
	double mean = calculateMean(size, noise);
	double std_dev = calculateStandardDeviation(size,noise);
	//printf("[");
	for (int i = 0; i < size; i++) {
		noise[i] = noise[i] - mean;
		noise[i] = noise[i] / (std_dev);      
		//printf("%lf, ", noise[i]);
	}
	//printf("]\n");

	free(fidx);
	free(faxis);
	free(filterf);
}