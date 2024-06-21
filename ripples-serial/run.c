#include "omislib.h"

#define ENOISE_BIN_PATH "C:\\Users\\mclab\\Desktop\\simone\\paral\\ripples-serial\\ripples-serial\\bin\\Enoise.bin"
#define INOISE_BIN_PATH "C:\\Users\\mclab\\Desktop\\simone\\paral\\ripples-serial\\ripples-serial\\bin\\Inoise.bin"


/*
*   @see https://en.wikipedia.org/wiki/Normal_distribution#Generating_values_from_normal_distribution
*/
double RandNormal()
{
    double x = 0, y = 0;
    while (x == 0 || y == 0) {
    x = (double)rand() / (double)RAND_MAX;
    y = (double)rand() / (double)RAND_MAX;
    }
    double z = sqrt(-2 * log(x)) * cos(2 * PI * y);

    return z;
}

/*
*   compare function (descending order) for quick sort
*/
static int DescendCmp(const void* a, const void* b) {
    double c = (*(double*)b - *(double*)a);
    if (c > 0)
        return -1;
    else if (c < 0)
        return 1;
    else
        return 0;
}


/*
* Calculate the matrix multiplication so that res = m1 * m2.
*/
void mulMatr(struct matrix* m1, struct matrix* m2, struct matrix *res) {
    long double sum;
    for (int i = 0; i < m1->size[0]; i++) {
        for (int j = 0; j < m2->size[1]; j++) {
            sum = 0;
            for (int k = 0; k < m1->size[1]; k++) {
                sum += m1->val[i * m1->size[1] + k] * m2->val[k * m2->size[1] + j];
            }
            res->val[i * res->size[1] + j] = (double)sum;
        }
    }
}

/*
* Calculate the matrix vector multiplication so that res = m * a.
* 
* @note: the number of rows of m must be equal to the size of a.
*/
void mulMatrVec(struct matrix* m, struct vec* a, double* res) {
    for (int i = 0; i < m->size[0]; i++) {
        res[i] = 0;
        for (int j = 0; j < m->size[1]; j++) {
            res[i] += m->val[i * m->size[1] + j] * a->val[j];
        }
    }
}

/*
* Calculate the product of two double arrays.
* 
* @param size: size of the arrays ( a -> [1 x size], b -> [size x 1] )
* @return a * b;
*/
double mulDArr(double *a, double *b, int size) {
    double sum = 0;
    for (int i = 0; i < size; i++) {
        sum += a[i] * b[i];
    }
    return sum;
}

/*
*   Calculate the transposed of a matrix.
*
*   @param mat The matrix.
*/
void transpose(struct matrix* mat) {
    struct matrix res;
    res.size[0] = mat->size[1];
    res.size[1] = mat->size[0];
    res.val = (double*)malloc(res.size[0] * res.size[1] * sizeof(double));
    for (int i = 0; i < mat->size[0]; i++) {
        for (int j = 0; j < mat->size[1]; j++) { 
            res.val[j * mat->size[0] + i] = mat->val[i * mat->size[1] + j];
        }
    }
    *mat = res;
}

void NetworkRunSeqt(struct pm p, struct inpseq in, int NE, int NI, double T, struct options opt) {

    /*outputs*/
    struct Conn conn;
    struct Vbar vbar;
    struct Veg veg;
    struct vec lfp;
    struct Tsp tspE;
    struct Tsp tspI;
    struct Isynbar isynbar;
    struct Inp inp;
    double* seqs;
    
    int max_size_tsp = 20000;
    tspE.times.size = max_size_tsp; // max size
    tspE.times.val = malloc(tspE.times.size * sizeof(double));
    tspE.celln.size = max_size_tsp; // max size
    tspE.celln.val = malloc(tspE.celln.size * sizeof(double));
    int tspE_count = 0;

    tspI.times.size = max_size_tsp; // max size
    tspI.times.val = malloc(tspI.times.size * sizeof(double));
    tspI.celln.size = max_size_tsp; // max size
    tspI.celln.val = malloc(tspI.celln.size * sizeof(double));
    int tspI_count = 0;

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

    free(w);


    /* inputs */
    struct matrix MX; //matrix to sum only some incoming inputs of CA3
    MX.size[0] = 100;
    MX.size[1] = NE;
    MX.val = (double*)malloc(MX.size[0] * MX.size[1] * sizeof(double));
    for (int i = 0; i < MX.size[0]; i++) {
        for (int j = 0; j < MX.size[1]; j++) {
            MX.val[i * MX.size[1] + j] = 0;
        }
    }
    int kkS;
    for (int i = 0; i < nL; i++) {
        kkS = floor(((float)NEseq[i] / NE) * 100) - 1;
        MX.val[kkS * MX.size[1] + NEseq[i]] = 1;
    }
    transpose(&MX);  // after it will be used transposed
    free(NEseq);

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
    transpose(&GEE); // after is used only transposed

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
    transpose(&GIE); // after is used only transposed

    clock_t toc = clock();
    printf("elapsed time is %.2lf seconds.\n", (double)(toc - tic) / CLOCKS_PER_SEC);


    /* initialize sim */


    //real dt is 0.001











    // time
    double dt = 0.001; // [=]ms integration step

    // allocate simulation ouput(GIE' * sIE)' * (vE - VrevI)
    int s = (ceil(1000 / dt) - 1) / 1000 + 1 + 1000 * (T - 1) + 1;
    vbar.E = (double*)malloc(s * sizeof(double));
    vbar.I = (double*)malloc(s * sizeof(double));
    veg.E = (double*)malloc(s * sizeof(double));
    veg.I = (double*)malloc(s * sizeof(double));
    lfp.size = s;
    lfp.val = (double*)malloc(s * sizeof(double));
    for (int i = 0; i < s; i++) {
        vbar.E[i] = 0;
        vbar.I[i] = 0;
        veg.E[i] = 0;
        veg.I[i] = 0;
        lfp.val[i] = 0;
    }

    // t = 0:dt:1000
    struct vec t; // one sec time axis 
    s = ceil(1000 / dt);
    t.size = s;
    t.val = (double*)malloc(s * sizeof(double));
    for (int i = 0; dt * i <= 1000; i++) {
        t.val[i] = dt * i;
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

    s = (ceil(1000 / dt) - 1) / 1000 + 1 + 1000 * (T - 1) + 1;
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

    s = ceil(1 / dt) + 1;
    struct matrix Einptrace;
    Einptrace.size[0] = NE;
    Einptrace.size[1] = s;
    Einptrace.val = (double*)malloc(Einptrace.size[0] * Einptrace.size[1] * sizeof(double));
    for (int i = 0; i < Einptrace.size[0]; i++) {
        for (int j = 0; j < Einptrace.size[1]; j++) {
            Einptrace.val[i * Einptrace.size[1] + j] = 0;
        }
    }

    s = ceil(1 / dt) + 1;
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
    Enoise.size[1] = (int)ceil(1000 / dt) + 1;
    //Enoise.size[1] = 301;
    struct matrix Inoise;
    Inoise.size[0] = NI;
    Inoise.size[1] = (int)ceil(1000 / dt) + 1;
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
    struct vec vE, wE, sE, sEI, erE, edE, erEI, edEI;
    vE.size = NE;
    wE.size = NE;
    sE.size = NE;
    sEI.size = NE;
    erE.size = NE;
    edE.size = NE;
    erEI.size = NE;
    edEI.size = NE;

    vE.val = (double *)malloc(NE * sizeof(double));
    wE.val = (double *)malloc(NE * sizeof(double));
    sE.val = (double *)malloc(NE * sizeof(double));
    sEI.val = (double *)malloc(NE * sizeof(double));
    erE.val = (double *)malloc(NE * sizeof(double));
    edE.val = (double *)malloc(NE * sizeof(double));
    erEI.val = (double *)malloc(NE * sizeof(double));
    edEI.val = (double *)malloc(NE * sizeof(double));

    struct vec vI, wI, sI, sIE, erI, edI, erIE, edIE;
    vI.size = NI;
    wI.size = NI;
    sI.size = NI;
    sIE.size = NI;
    erI.size = NI;
    edI.size = NI;
    erIE.size = NI;
    edIE.size = NI;

    vI.val = (double *)malloc(NI * sizeof(double));
    wI.val = (double *)malloc(NI * sizeof(double));
    sI.val = (double *)malloc(NI * sizeof(double));
    sIE.val = (double *)malloc(NI * sizeof(double));
    erI.val = (double *)malloc(NI * sizeof(double));
    edI.val = (double *)malloc(NI * sizeof(double));
    erIE.val = (double *)malloc(NI * sizeof(double));
    edIE.val = (double *)malloc(NI * sizeof(double));


    double tmin, tmax;
    struct vec stsec;
    stsec.size = in.on.size; //max size, then realloc
    stsec.val = (double *)malloc(stsec.size * sizeof(double));

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
        if (seqN == 1) {
            for (int i = 0; i < NE; i++) {
                vE.val[i] = ((double)rand() / (double)RAND_MAX) * (70 + p.VrE) - 70;
                wE.val[i] = p.aE * (vE.val[i] - p.ElE);
                sE.val[i] = 0;
                sEI.val[i] = 0;
                erE.val[i] = 0;
                edE.val[i] = 0;
                erEI.val[i] = 0;
                edEI.val[i] = 0;
            }

            for (int i = 0; i < NI; i++) {
                vI.val[i] = ((double)rand() / (double)RAND_MAX) * (70 + p.VrI) - 70;
                wI.val[i] = p.aI * (vI.val[i] - p.ElI);
                sI.val[i] = 0;
                sIE.val[i] = 0;
                erI.val[i] = 0;
                edI.val[i] = 0;
                erIE.val[i] = 0;
                edIE.val[i] = 0;
            }
        }
        else {
            tmin = (seqN - 1) * 1000 - 100;
            tmax = seqN * 1000 + 20;
            stsec.size = 0;
            for (int i = 0; i < in.on.size; i++) {
                if (in.on.val[i] >= tmin && in.on.val[i] < tmax) {
                    stsec.val[stsec.size] = in.on.val[i]; // note that in the matlab version this part is wrong
                    stsec.size++;
                }
            }
            stsec.val = realloc(stsec.val, stsec.size * sizeof(double));

            struct matrix bmps;
            double *bt, *ebt;
            double ebt_max = 0;
            double bt_max = 0;
            bmps.size[0] = 100;
            bmps.size[1] = t.size;
            bmps.val = (double *)malloc(bmps.size[0] * bmps.size[1] * sizeof(double));
            bt = (double*)malloc(t.size * sizeof(double));
            ebt = (double*)malloc(t.size * sizeof(double));

            for (int j = 0; j < t.size; j++) {
                for (int i = 0; i < 100; i++) {
                    bmps.val[i * bmps.size[1] + j] = 0;
                }
                bt[j] = 0;
                ebt[j] = 0;
            }

            for (int i = 0; i < stsec.size; i++) {
                // inside each ripple event
                stsec.val[i] = stsec.val[i] - ((seqN - 1) * 1000);
                double rplton = stsec.val[i];
                double rpltoff = rplton + in.length;
                int L = in.length - in.slp - 2;
                int L0 = 0; // in.slp/in.length
                int L1 = 1; // - L0

                double step = (double)1 / 99;
                int size = ceil((L1 - L0) / step) + 1;
                double* tbins, *tons, *toffs;
                tbins = (double*)malloc(size * sizeof(double));
                tons = (double*)malloc(size * sizeof(double));
                toffs = (double*)malloc(size * sizeof(double));
                for (int k = 0; k < size; k++) {
                    tbins[k] = (L0 + (k * step)) * L;
                    tons[k] = rplton + in.slp + 2 + tbins[k]; // start the Ecells bumps after the I cells are inhibiting already
                    toffs[k] = tons[k] + (in.length / 99);
                }
                for (int j = 0; j < bmps.size[1]; j++) {
                    for (int k = 0; k < bmps.size[0]; k++) {
                        bmps.val[k * bmps.size[1] + j] = bmps.val[k * bmps.size[1] + j] + ( 1 / (1 + exp(( tons[k] - t.val[j]) / 1.5)) * 1 / (1 + exp((t.val[j] - toffs[k]) / 1.5)));
                    }
                    bt[j] = bt[j] + ( (1 / (1 + exp((rplton - t.val[j]) / in.slp))) * (1 / (1 + exp((t.val[j] - rpltoff) / in.slp))));
                    ebt[j] = ebt[j] + ( (1 / (1 + exp((rplton - t.val[j]) / (in.slp/2)))) * (1 / (1 + exp((t.val[j] - rpltoff) / (in.slp/2)))));
                    ebt_max = ebt_max > ebt[j] ? ebt_max : ebt[j];
                    bt_max = bt_max > bt[j] ? bt_max : bt[j];
                }

                free(tbins);
                free(tons);
                free(toffs);
            }
            
            struct matrix AEX, AIX;
            AEX.size[0] = NE;
            AEX.size[1] = bmps.size[1];
            AIX.size[0] = NI;
            AIX.size[1] = bmps.size[1];
            AEX.val = (double*)malloc(AEX.size[0] * AEX.size[1] * sizeof(double));
            AIX.val = (double*)malloc(AIX.size[0] * AIX.size[1] * sizeof(double));

            struct matrix tmp;
            tmp.size[0] = MX.size[0];
            tmp.size[1] = bmps.size[1];
            tmp.val = (double*)malloc(tmp.size[0] * tmp.size[1] * sizeof(double));
            mulMatr(&MX, &bmps, &tmp);
            
            for (int k = 0; k < AEX.size[0]; k++) {
                for (int j = 0; j < AEX.size[1]; j++) {
                    AEX.val[k * AEX.size[1] + j] = 5 * tmp.val[k * tmp.size[1] + j] + ebt[j];
                }
            }

            for (int k = 0; k < AIX.size[0]; k++) {
                for (int j = 0; j < AIX.size[1]; j++) {
                    AIX.val[k * AIX.size[1] + j] = bt[j];
                }
            }
            free(tmp.val);

            double Escale = ebt_max;
            double ff;
            ff = Escale > 0 ? p.jmpE / Escale : 0;
            for (int k = 0; k < Enoise.size[0]; k++) {
                for (int j = 0; j < Enoise.size[1]; j++) {
                    Enoise.val[k * Enoise.size[1] + j] = Enoise.val[k * Enoise.size[1] + j] + (ff * AEX.val[k * AEX.size[1] + j]);
                }
            }            
            double Iscale = bt_max;
            double gg;
            gg = Iscale > 0 ? p.jmpI / Iscale : 0;
            for (int k = 0; k < Inoise.size[0]; k++) {
                for (int j = 0; j < Inoise.size[1]; j++) {
                    Inoise.val[k * Inoise.size[1] + j] = Inoise.val[k * Inoise.size[1] + j] + (gg * AIX.val[k * Inoise.size[1] + j]);
                }
            }
            
            int step = 1 / dt;
            Einptrace.size[1] += AEX.size[1] / step; // should be always 1000 (so += AEX.size[1] / step) or lower if dt is bigger? (so += step)
            Einptrace.val = realloc(Einptrace.val, Einptrace.size[0] * Einptrace.size[1] * sizeof(double));
            for (int k = 0; k < Einptrace.size[0]; k++) {
                for (int j = step+1; j < Einptrace.size[1]; j++) {
                    Einptrace.val[k * Einptrace.size[1] + j] =  ff * AEX.val[k * AEX.size[1] + ((j - step) * step)];
                }
            }            
            Iinptrace.size[1] += AIX.size[1] / step; // should be always 1000 (so += AEX.size[1] / step) or lower if dt is bigger? (so += step)
            Iinptrace.val = realloc(Iinptrace.val, Iinptrace.size[0] * Iinptrace.size[1] * sizeof(double));
            for (int k = 0; k < Iinptrace.size[0]; k++) {
                for (int j = step+1; j < Iinptrace.size[1]; j++) {
                    Iinptrace.val[k * Iinptrace.size[1] + j] =  gg * AIX.val[k * AIX.size[1] + ((j - step) * step)];
                }
            }


            free(bmps.val);
            free(bt);
            free(ebt);
            free(AEX.val); 
            
        }
        
        int ids;

        double* tmp1, *tmp2, * tmp3, * tmp4;
        tmp1 = (double*)malloc(NE * sizeof(double));
        tmp2 = (double*)malloc(NE * sizeof(double));
        tmp3 = (double*)malloc(NI * sizeof(double));
        tmp4 = (double*)malloc(NE * sizeof(double));

        mulMatrVec(&GEE, &sE, tmp1); 
        mulMatrVec(&GIE, &sIE, tmp2); 
        mulMatrVec(&GII, &sI, tmp3); 
        mulMatrVec(&GEI, &sEI, tmp4); 


        double Isyn, Iapp, Iion, dvdt, dwdt, wEst, wIst, vEst, vIst;

        for (int idt = 1; idt < t.size; idt++) {
            if (idt % 1000 == 1 && idt > 1) { 
                ids = (idt - 1) / 1000 + 1 + 1000 * (seqN - 1);
                double vE_sum = 0, vI_sum = 0;
                for (int i = 0; i < NE; i++) {
                    if (vE.val[i] > 0)
                        vE.val[i] = 50;
                    vE_sum += vE.val[i];
                }
                for (int i = 0; i < NI; i++) {
                    if (vI.val[i] > 0)
                        vI.val[i] = 50;
                    vI_sum += vI.val[i];
                }
                vbar.E[ids] = vE_sum / NE;
                vbar.I[ids] = vI_sum / NI;
                veg.E[ids] = vE.val[Eeg];
                veg.I[ids] = vI.val[Ieg];
                lfp.val[ids] = 0;
                for (int i = 0; i < NE; i++) {
                    lfp.val[ids] += ( tmp1[i] * (vE.val[i] - p.VrevE) ) + (tmp2[i] * (vE.val[i] - p.VrevI)); 
                }

                lfp.val[ids] /= NE;

                if (opt.storecurrs) {
                    for (int i = 0; i < NE; i++) {
                        isynbar.EtoE.val[i * isynbar.EtoE.size[1] + ids] = tmp1[i] * (vE.val[i] - p.VrevE);
                        isynbar.ItoE.val[i * isynbar.ItoE.size[1] + ids] = tmp2[i] * (vE.val[i] - p.VrevI);
                    }
                }
                
            }

            /* E cells */
            for (int idx = 0; idx < NE; idx++) {
                vEst = vE.val[idx];
                wEst = wE.val[idx];
                Isyn = (tmp1[idx] * (vE.val[idx] - p.VrevE)) + (tmp2[idx] * (vE.val[idx] - p.VrevI) );
                Iapp = Enoise.val[idx * Enoise.size[1] + idt];
                Iion = ((-1) * p.glE * (vE.val[idx] - p.ElE)) + (p.glE * p.slpE * exp((vE.val[idx] - p.VtE) / p.slpE)) - (wE.val[idx]);
                
                dvdt = ((Iapp + Iion - Isyn) / p.CE);
                dwdt = ((p.aE * (vE.val[idx] - p.ElE) - wE.val[idx]) / p.twE);
                vE.val[idx] = vE.val[idx] + (dt * dvdt);
                wE.val[idx] = wE.val[idx] + (dt * dwdt);

                // syn gates evolution
                edE.val[idx] *= fdE;  
                erE.val[idx] *= frE;
                sE.val[idx] = erE.val[idx] - edE.val[idx];
                edEI.val[idx] *= fdEI;
                erEI.val[idx] *= frEI;
                sEI.val[idx] = erEI.val[idx] - edEI.val[idx];
                
                
                if (vEst >= 0) {
                    // update dynamic vars
                    vE.val[idx] = p.VrE;
                    wE.val[idx] = wEst + p.bE;
                    // update syn gates
                    edE.val[idx] += 1 / pvsE;
                    erE.val[idx] += 1 / pvsE;
                    edEI.val[idx] += 1 / pvsEI;
                    erEI.val[idx] += 1 / pvsEI;

                    tspE.times.val[tspE_count] = (t.val[idt] / 1000) + seqN - 1;
                    tspE.celln.val[tspE_count] = idx;
                    tspE_count++;
                }
            }

            /* I cells */
            for (int idx = 0; idx < NI; idx++) {
                vIst = vI.val[idx];
                wIst = wI.val[idx];
                Isyn = (tmp3[idx] * (vI.val[idx] - p.VrevI)) + (tmp4[idx] * (vI.val[idx] - p.VrevE));
                Iapp = Inoise.val[idx * Enoise.size[1] + idt];
                Iion = ((-1) * p.glI * (vI.val[idx] - p.ElI)) + (p.glI * p.slpI * exp((vI.val[idx] - p.VtI) / p.slpI)) - (wI.val[idx]);

                dvdt = ((Iapp + Iion - Isyn) / p.CI);
                dwdt = ((p.aI * (vI.val[idx] - p.ElI) - wI.val[idx]) / p.twI);
                vI.val[idx] += dt * dvdt;
                wI.val[idx] += dt * dwdt;

                // syn gates evolution
                edI.val[idx] *= fdI;
                erI.val[idx] *= frI;
                sI.val[idx] = erI.val[idx] - edI.val[idx];
                edIE.val[idx] *= fdIE;
                erIE.val[idx] *= frIE;
                sIE.val[idx] = erIE.val[idx] - edIE.val[idx];


                if (vIst >= 0) { 
                    // update dynamic vars
                    vI.val[idx] = p.VrI;
                    wI.val[idx] = wIst + p.bI;
                    // update syn gates
                    edI.val[idx] += 1 / pvsI;
                    erI.val[idx] += 1 / pvsI;
                    edIE.val[idx] += 1 / pvsIE;
                    erIE.val[idx] += 1 / pvsIE;

                    tspI.times.val[tspI_count] = (t.val[idt] / 1000) + seqN - 1;
                    tspI.celln.val[tspI_count] = idx;
                    tspI_count++;
                }
            }

        }

        free(tmp1);
        free(tmp2);

        seqN++;
        clock_t toc = clock();
        printf("elapsed time is %.2lf seconds.\n", (double)(toc - tic) / CLOCKS_PER_SEC);
    }
    
    tspE.times.size = tspE_count;
    tspE.celln.size = tspE_count;
    tspE.times.val = realloc(tspE.times.val, tspE.times.size * sizeof(double));
    tspE.celln.val = realloc(tspE.celln.val, tspE.celln.size * sizeof(double));   
    
    tspI.times.size = tspI_count;
    tspI.celln.size = tspI_count;
    tspI.times.val = realloc(tspI.times.val, tspI.times.size * sizeof(double));
    tspI.celln.val = realloc(tspI.celln.val, tspI.celln.size * sizeof(double));

    inp.Etrace = Einptrace;
    inp.Itrace = Iinptrace;


    
    /* free */
    free(MX.val);
    free(GEE.val);
    free(GII.val);
    
    free(Enoise.val);
    free(Inoise.val);

    free(vE.val); 
    free(wE.val); 
    free(sE.val); 
    free(sEI.val); 
    free(erE.val); 
    free(edE.val); 
    free(erEI.val); 
    free(edEI.val);
    free(vI.val);
    free(wI.val);
    free(sI.val);
    free(sIE.val);
    free(erI.val);
    free(edI.val);
    free(erIE.val);
    free(edIE.val);
    free(stsec.val);
    
    return;
}