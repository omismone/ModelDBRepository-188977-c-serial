#include "omislib.h"

#define ENOISE_BIN_PATH "C:\\Users\\mclab\\Desktop\\simone\\paral\\ripples-serial\\ripples-serial\\bin\\Enoise.bin"
#define INOISE_BIN_PATH "C:\\Users\\mclab\\Desktop\\simone\\paral\\ripples-serial\\ripples-serial\\bin\\Inoise.bin"


/*
*   @see https://en.wikipedia.org/wiki/Normal_distribution#Generating_values_from_normal_distribution
*/
double RandNormal()
{
    double x = (double)rand() / (double)RAND_MAX;
    double y = (double)rand() / (double)RAND_MAX;
    double z = sqrt(-2 * log(x)) * cos(2 * PI * y);

    return z;
}

/*
*   compare function (descending order) for quick sort
*/
int DescendCmp(const void* a, const void* b) {
    double c = (*(double*)b - *(double*)a);
    if (c > 0)
        return -1;
    else if (c < 0)
        return 1;
    else
        return 0;
}


void NetworkRunSeqt(struct pm p, struct inpseq in, int NE, int NI, float T, struct options opt) {

    /*outputs*/
    struct Conn conn;
    struct Vbar vbar;
    struct Veg veg;
    double* lfp;
    struct Tsp tspE;
    struct Tsp tspI;
    struct Isynbar isynbar;
    struct Inp inp;
    double* seqs;

    /* pre select the sequence */
    int* NEseq;
    int nL = 10;
    NEseq = (int*)malloc(nL * sizeof(int));
    for (int i = 0; i < nL; i++) {
        NEseq[i] = ((rand() % 80) + 1) + (i * 80);
    }



    /* opt options: nonoise novar noiseprc */
    p.gnoiseE = opt.nonoise ? 0 : p.gnoiseE * (opt.noiseprc / 100);
    p.gnoiseI = opt.nonoise ? 0 : p.gnoiseI * (opt.noiseprc / 100);

    double* Edc_dist;
    Edc_dist = (double*)malloc(NE * sizeof(double));
    double Edc_dist_max = log(0); //= -inf
    for (int i = 0; i < NE; i++) {
        Edc_dist[i] = (RandNormal() * p.DCstdE * p.Edc) + p.Edc;
        Edc_dist_max = Edc_dist[i] > Edc_dist_max ? Edc_dist[i] : Edc_dist_max;
    }
    double* Idc_dist;
    Idc_dist = (double*)malloc(NI * sizeof(double));
    for (int i = 0; i < NI; i++) {
        Idc_dist[i] = (RandNormal() * p.DCstdI * p.Idc) + p.Idc;
    }

    double* w;
    w = (double*)malloc(nL * sizeof(double));
    for (int i = 0; i < nL; i++) {
        w[i] = RandNormal();
    }
    qsort(w, nL, sizeof(double), DescendCmp);
    for (int i = 0; i < nL; i++) {
        Edc_dist[NEseq[i]] = (Edc_dist_max + (p.dcbias * p.DCstdE * p.Edc)) * ((double)1 + (w[i] * (double)0.02));
    }
    inp.Edc = Edc_dist; //on matlab is transposed
    inp.Idc = Idc_dist;



    /* inputs */
    struct matrix MX; //matrix to sum only some incoming inputs of CA3
    MX.size[0] = 100;
    MX.size[1] = NE;
    MX.val = (double*)malloc(MX.size[0] * MX.size[1] * sizeof(double));
    int kkS;
    for (int i = 0; i < nL; i++) {
        kkS = ceil((NEseq[i] / NE) * 100);
        MX.val[kkS * MX.size[1] + NEseq[i]] = 1;
    }

    /* synapses */
    printf("wire ntwk - ");
    clock_t tic = clock();
    if (opt.novar) {
        p.gvarEE = 0;
        p.gvarII = 0;
        p.gvarEI = 0;
        p.gvarIE = 0;
    }

    double mn = p.gmaxEE / NE;
    double vr = p.gvarEE * mn;
    struct matrix GEE;
    GEE.size[0] = NE;
    GEE.size[1] = NE;
    GEE.val = (double*)malloc(GEE.size[0] * GEE.size[1] * sizeof(double));
    double g;
    for (int i = 0; i < GEE.size[0]; i++) {
        for (int j = 0; j < GEE.size[1]; j++) {
            g = RandNormal() * sqrt(vr) + mn;
            GEE.val[i * GEE.size[1] + j] = g > 0 ? g : 0;
        }
    }
    conn.EtoE = GEE;

    mn = p.gmaxII / NI;
    vr = p.gvarII * mn;
    struct matrix GII;
    GII.size[0] = NI;
    GII.size[1] = NI;
    GII.val = (double*)malloc(GII.size[0] * GII.size[1] * sizeof(double));
    for (int i = 0; i < GII.size[0]; i++) {
        for (int j = 0; j < GII.size[1]; j++) {
            g = RandNormal() * sqrt(vr) + mn;
            GII.val[i * GII.size[1] + j] = g > 0 ? g : 0;
        }
    }
    conn.ItoI = GII;

    mn = p.gmaxEI / NE;
    vr = p.gvarEI * mn;
    struct matrix GEI;
    GEI.size[0] = NE;
    GEI.size[1] = NI;
    GEI.val = (double*)malloc(GEI.size[0] * GEI.size[1] * sizeof(double));
    for (int i = 0; i < GEI.size[0]; i++) {
        for (int j = 0; j < GEI.size[1]; j++) {
            g = RandNormal() * sqrt(vr) + mn;
            GEI.val[i * GEI.size[1] + j] = g > 0 ? g : 0;
        }
    }
    conn.EtoI = GEI;

    mn = p.gmaxIE / NI;
    vr = p.gvarIE * mn;
    struct matrix GIE;
    GIE.size[0] = NI;
    GIE.size[1] = NE;
    GIE.val = (double*)malloc(GIE.size[0] * GIE.size[1] * sizeof(double));
    for (int i = 0; i < GIE.size[0]; i++) {
        for (int j = 0; j < GIE.size[1]; j++) {
            g = RandNormal() * sqrt(vr) + mn;
            GIE.val[i * GIE.size[1] + j] = g > 0 ? g : 0;
        }
    }
    conn.ItoE = GIE;

    clock_t toc = clock();
    printf("elapsed time is %.2lf seconds.\n", (double)(toc - tic) / CLOCKS_PER_SEC);


    /* initialize sim */

    // allocate simulation ouput(GIE' * sIE)' * (vE - VrevI)
    int s = ceil(T / 0.001);
    vbar.E = (double*)malloc(s * sizeof(double));
    vbar.I = (double*)malloc(s * sizeof(double));
    veg.E = (double*)malloc(s * sizeof(double));
    veg.I = (double*)malloc(s * sizeof(double));
    lfp = (double*)malloc(s * sizeof(double));
    for (int i = 0; i < s; i++) {
        vbar.E[i] = 0;
        vbar.I[i] = 0;
        veg.E[i] = 0;
        veg.I[i] = 0;
        lfp[i] = 0;
    }

    // time
    double dt = 0.001; // [=]ms integration step
    // t = 0:dt:1000
    double* t; // one sec time axis 
    s = ceil(1000 / dt);
    t = (double*)malloc(s * sizeof(double));
    for (int i = 0; dt * i <= 1000; i++) {
        t[i] = dt * i;
    }

    // peak values of biexps signals
    double pvsE = exp((1 / (1 / p.tauEd - 1 / p.tauEr) * log(p.tauEd / p.tauEr)) / p.tauEr) - exp((1 / (1 / p.tauEd - 1 / p.tauEr) * log(p.tauEd / p.tauEr)) / p.tauEd);
    double fdE = exp(-dt / p.tauEd); //factor of decay
    double frE = exp(-dt / p.tauEr); //factor of rise

    double pvsI = exp((1 / (1 / p.tauId - 1 / p.tauIr) * log(p.tauId / p.tauIr)) / p.tauIr) - exp((1 / (1 / p.tauId - 1 / p.tauIr) * log(p.tauId / p.tauIr)) / p.tauId);
    double fdI = exp(-dt / p.tauId);
    double frI = exp(-dt / p.tauIr);

    double pvsIE = exp((1 / (1 / p.tauIEd - 1 / p.tauIEr) * log(p.tauIEd / p.tauIEr)) / p.tauIEr) - exp((1 / (1 / p.tauIEd - 1 / p.tauIEr) * log(p.tauIEd / p.tauIEr)) / p.tauIEd);
    double fdIE = exp(-dt / p.tauIEd);
    double frIE = exp(-dt / p.tauIEr);

    double pvsEI = exp((1 / (1 / p.tauEId - 1 / p.tauEIr) * log(p.tauEId / p.tauEIr)) / p.tauEIr) - exp((1 / (1 / p.tauEId - 1 / p.tauEIr) * log(p.tauEId / p.tauEIr)) / p.tauEId);
    double fdEI = exp(-dt / p.tauEId);
    double frEI = exp(-dt / p.tauEIr);

    int Eeg = (rand() % NE) + 1;
    veg.ne = Eeg;
    int Ieg = (rand() % NI) + 1;
    veg.ni = Ieg;

    s = ceil(T / 0.001) + 1;
    isynbar.ItoE.size[0] = NE;
    isynbar.ItoE.size[1] = s;
    isynbar.ItoE.val = (double*)malloc(isynbar.ItoE.size[0] * isynbar.ItoE.size[1] * sizeof(double));
    isynbar.EtoE.size[0] = NE;
    isynbar.EtoE.size[1] = s;
    isynbar.EtoE.val = (double*)malloc(isynbar.EtoE.size[0] * isynbar.EtoE.size[1] * sizeof(double));
    if (opt.storecurrs) {
        for (int i = 0; i < isynbar.ItoE.size[0]; i++) {
            for (int j = 0; j < isynbar.ItoE.size[1]; j++) {
                isynbar.ItoE.val[i * isynbar.ItoE.size[1] + j] = 0;
            }
        }
        for (int i = 0; i < isynbar.EtoE.size[0]; i++) {
            for (int j = 0; j < isynbar.EtoE.size[1]; j++) {
                isynbar.EtoE.val[i * isynbar.EtoE.size[1] + j] = 0;
            }
        }
    }

    s = ceil(1 / 0.001) + 1;
    struct matrix Einptrace;
    Einptrace.size[0] = NE;
    Einptrace.size[1] = s;
    Einptrace.val = (double*)malloc(Einptrace.size[0] * Einptrace.size[1] * sizeof(double));
    for (int i = 0; i < Einptrace.size[0]; i++) {
        for (int j = 0; j < Einptrace.size[1]; j++) {
            Einptrace.val[i * Einptrace.size[1] + j] = 0;
        }
    }

    s = ceil(1 / 0.001) + 1;
    struct matrix Iinptrace;
    Iinptrace.size[0] = NI;
    Iinptrace.size[1] = s;
    Iinptrace.val = (double*)malloc(Iinptrace.size[0] * Iinptrace.size[1] * sizeof(double));
    for (int i = 0; i < Iinptrace.size[0]; i++) {
        for (int j = 0; j < Iinptrace.size[1]; j++) {
            Iinptrace.val[i * Iinptrace.size[1] + j] = 0;
        }
    }

    /* for each second of simulation we generate new noise */

    /* the noise is taken by a binary file, alternatively it can be generated with the function noiseGen but the size has to be small and odd (e.g: 301) */
   
    int seqN = 1;
    struct matrix Enoise;
    Enoise.size[0] = NE;
    Enoise.size[1] = (int)ceil(1000 / 0.001) + 1;
    //Enoise.size[1] = 301;
    struct matrix Inoise;
    Inoise.size[0] = NI;
    Inoise.size[1] = (int)ceil(1000 / 0.001) + 1;
    //Inoise.size[1] = 301;
    Enoise.val = (double*)malloc(Enoise.size[0] * Enoise.size[1] * sizeof(double));
    Inoise.val = (double*)malloc(Inoise.size[0] * Inoise.size[1] * sizeof(double));

    // [CHOISE 1 - recommended] reading the noise, to have it in a binary file @see scripts/SaveNoise.m 
    printf("[reading noise]\n");

    FILE* fp = NULL;
    fp = fopen(ENOISE_BIN_PATH, "rb");
    if (fp == NULL) {
        perror("error while opening noises files");
        exit(1);
    }
    int numElements = fread(Enoise.val, sizeof(double), Enoise.size[0] * Enoise.size[1], fp);
    if (numElements < Enoise.size[0] * Enoise.size[1]) {
        if (ferror(fp)) {
            perror("error while reading noises files");
        }
        else {
            printf("[[warning: EOF before reading all the noises]]\n");
        }
        fclose(fp);
        exit(1);
    }
    fclose(fp);

    fp = fopen(INOISE_BIN_PATH, "rb");
    if (fp == NULL) {
        perror("error while opening noises files");
        exit(1);
    }
    numElements = fread(Inoise.val, sizeof(double), Inoise.size[0] * Inoise.size[1], fp);
    if (numElements < Inoise.size[0] * Inoise.size[1]) {
        if (ferror(fp)) {
            perror("error while reading noises files");
        }
        else {
            printf("[[warning: EOF before reading all the noises]]\n");
        }
        fclose(fp);
        exit(1);
    }
    fclose(fp);


    // for the initialization of the first instant of simulation
    double* vE, * wE, * sE, * sEI, * erE, * edE, * erEI, * edEI;
    vE = malloc(NE * sizeof(double));
    wE = malloc(NE * sizeof(double));
    sE = malloc(NE * sizeof(double));
    sEI = malloc(NE * sizeof(double));
    erE = malloc(NE * sizeof(double));
    edE = malloc(NE * sizeof(double));
    erEI = malloc(NE * sizeof(double));
    edEI = malloc(NE * sizeof(double));

    double* vI, * wI, * sI, * sII, * erI, * edI, * erIE, * edIE;
    vI = malloc(NI * sizeof(double));
    wI = malloc(NI * sizeof(double));
    sI = malloc(NI * sizeof(double));
    sII = malloc(NI * sizeof(double));
    erI = malloc(NE * sizeof(double));
    edI = malloc(NE * sizeof(double));
    erIE = malloc(NE * sizeof(double));
    edIE = malloc(NE * sizeof(double));

    while (seqN <= T) {
        /*
        // [CHOISE 2] generating the noise 

        tic = clock();
        printf("[noise generation] second %d\n", seqN);
        for (int i = 0; i < Enoise.size[0]; i++) { //the noise is generated row by row
            noiseGen(&Enoise.val[i* Enoise.size[1]], Enoise.size[1], 1, dt);
        }
        for (int i = 0; i < Inoise.size[0]; i++) {
            noiseGen(&Inoise.val[i* Inoise.size[1]], Inoise.size[1], 1, dt);
        }
        */

        printf("integrating ODE\n");


        for (int i = 0; i < NE; i++) {
            vE[i] = ((double)rand() / (double)RAND_MAX) * (70+p.VrE) - 70;
            wE[i] = p.aE * (vE[i] - p.ElE);
            sE[i] = 0;
            sEI[i] = 0;
            if (seqN == 1) {
                erE[i] = 0;
                edE[i] = 0;
                erEI[i] = 0;
                edEI[i] = 0;
            }
        }

        for (int i = 0; i < NI; i++) {
            vI[i] = ((double)rand() / (double)RAND_MAX) * (70 + p.VrI) - 70;
            wI[i] = p.aI * (vI[i] - p.ElI);
            sI[i] = 0;
            sII[i] = 0;
            if (seqN == 1) {
                erI[i] = 0;
                edI[i] = 0;
                erIE[i] = 0;
                edIE[i] = 0;
            }
        }
        int tmin, tmax;
        if (seqN != 1) {
            tmin = (seqN - 1) * 1000 - 100;
            tmax = seqN * 1000 + 20;

            // to do: comnporre il vettore stsec
        }
        


        seqN++;
    }




    /* free */
    free(NEseq);
    free(Edc_dist);
    free(Idc_dist);
    free(w);
    free(MX.val);
    free(GEE.val);
    free(GII.val);
    // free(vbar.E);
    // free(vbar.I);
    // free(veg.E);
    // free(veg.I);
    // free(lfp);
    free(Enoise.val);
}