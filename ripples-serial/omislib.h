#pragma once
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include<complex.h> 
#undef I
#define J _Complex_I

#define PI 3.14159265358979323846

/*
*   Array structure that contains the array and its size.
*/
struct vec {
    double* val;   // The array.
    int size;    // Number of its elements.
};

/*
*   Array structure that contains the array and its size.
*/
struct long_vec {
    double* val;   // The array.
    int size;    // Number of its elements.
};

/*
*   Matrix structure that contains the matrix and its size. //optimizing memory access.
*/
struct matrix {
    double* val;   // The values.
    int size[2];    // size[0] # of rows, size[1] # of cols.
};

/*
*   Hippocampal pyramidal layer parameters
*/
struct pm {
    /* E = excitatory ; I = inhibitory */

    double CE;   // capacity [F]
    double glE;  // leakage conductance [S]
    double ElE;  // leakage reversal potential [V]
    double aE;   // constant [ ]
    double bE;   // constant [ ]
    double slpE; // 
    double twE;  // [s]
    double VtE;  // threshold
    double VrE;  // resting voltage

    double CI;
    double glI;
    double ElI;
    double aI;
    double bI;
    double slpI;
    double twI;
    double VtI;
    double VrI;

    double gnoiseE;
    double gnoiseI;

    //EtoE
    double tauEr;
    double tauEd;
    //EtoB
    double tauEIr;
    double tauEId;
    //BtoB
    double tauIr;
    double tauId;
    //BtoP 
    double tauIEr;
    double tauIEd;

    double gmaxII;
    double gmaxEI;
    double gmaxIE;
    double gmaxEE;

    double VrevE;
    double VrevI;

    double gvarEE;
    double gvarII;

    double gvarEI;
    double gvarIE;

    double DCstdI;
    double DCstdE;

    double Edc;
    double jmpE;
    double Idc;
    double jmpI;

    double seqsize;
    double dcbias;
};

/*
*   input sequence
*/
struct inpseq {
    float slp;
    struct vec on;
    float length;
};

struct options {
    int nonoise; //if no noise added, turn to 1
    int novar; //if no variance in synaptic weightsm turn to 1
    int storecurrs; //if you want the output to include the synaptic currents
    float noiseprc; //percent of standard deviation of the noise to use in the simulation
    int seqassign; //if you want to choose 10 cells that are going to be part of a sequence
};

/*
*   neurons connections
*/
struct Conn {
    struct matrix EtoE;
    struct matrix ItoI;
    struct matrix EtoI;
    struct matrix ItoE;
};

/*
*   Y_syn v_n
*/
struct Vbar {
    double* E;
    double* I;
};

/*
*   outputs
*/
struct Veg {
    double* E;
    double* I;
    int ne;
    int ni;
};

/*
*  spikes times
*/
struct Tsp {
    struct vec times;
    struct vec celln;
};

/*
*   I_syn
*/
struct Isynbar {
    struct matrix ItoE;
    struct matrix EtoE;
};

/*
*   don't need this if all neurons have the same inputs
*/
struct Inp {
    double* Edc;
    double* Idc;
    struct matrix Etrace;
    struct matrix Itrace;
};


/*
*   run the simulation.
*/
void NetworkRunSeqt(struct pm p, struct inpseq in, int NE, int NI, double T, struct options opt);

/*
*  generate the UO noise.
* 
* @param noise: where to save the generated noise vector.
* @param size: the size of the noise vector.
* @param T: the duration (?)
*/
void noiseGen(double* noise, int size, double T, double dt);

double RandNormal();

