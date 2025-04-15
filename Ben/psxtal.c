/*
 *  xtal.c
 *
 *
 *  Created by marisol on 11/17/08.
 *  Copyright 2008 Purdue. All rights reserved.
 *
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <dfftw_mpi.h>
#include <string.h>

#define DEFAULTV -200
short N1		= DEFAULTV;
short N2 	    = DEFAULTV;
short N3 		= DEFAULTV;
int   NT 		= DEFAULTV;
short load		= DEFAULTV;
short NSI 		= DEFAULTV;
short itp 	    = DEFAULTV;
short rsm		= DEFAULTV;
short NS		= DEFAULTV;   /* 12 number of slip systems for FCC */
short ND		= DEFAULTV;
short NP		= DEFAULTV;
short iflag		= DEFAULTV;
short nloops	= DEFAULTV;
short config	= DEFAULTV;
short radius2	= DEFAULTV;
short mflag	= DEFAULTV;
short cflag	= DEFAULTV;
float An	= DEFAULTV;
float Bn	= DEFAULTV;
short oflag	= DEFAULTV;
float setsigma	= DEFAULTV;
float sigstep	= DEFAULTV;
float seteps    = DEFAULTV;
float epstep    = DEFAULTV;
char scratch[100];
short voro	= DEFAULTV;
short test	= DEFAULTV;
char StandardFile[100];
short dump_3 = DEFAULTV;
short comment = DEFAULTV;


// 2mat Cu-Ni elasticity stress
#define CuX 3408650
#define CuY 72219.171875
#define NiX	5845149.5
#define NiY	-72219.679688

void 	inputVariables(char *PFDDinput);
void	checkVariables();

double **epsilon;

int main(int argc, char **argv)
{
    inputVariables("PFDD_input.dat"); // Call a function to input the variables from an extern file
    
    int i, j, k, m, n, ksym, it, tmstep, is, is2, id, ir, ka, kb, nb, na, ida, idb;
    int na0, na1, na2, nsize, k1, k2, k3, nfreq, ii, b1, b2;
    int isa, isb, psys, fflag;
    int rank, numprocs, rtn_val;
    int lN1, lxs, lnyt, lyst, lsize;
    int index, index1, index2, ia, ib, ic, ig, ie, ih, ij;
    int ic1, ic2, ic3, ic4;
    unsigned long seed;
    double *fx, *fy, *fz, *xi, *xi_sum, *xo, pi;
    double *f, *r, *fcore, *df1core, *df2core, *df3core;
    double L, d1, d2, d3, dt, size, size3, sizegs, xi_ave, xiave, alpha;
    double *dE_core, E_core;
    double *data_eps, *data_epsd, *data_sigma;
    double c0, c1, c2, c3, c4, a1, a3, isf, usf; /*constants for gamma surface*/
    double b, C11, C12, C44, dslip, CD, obsden, dtt, mu, young, nu, a, ll, C[3][3][3][3];
    double *BB, *FF, *DD, *UU, S11, S12, S44;
    double T0, TT0, T1, TT1, T2, TT2, T9, TT9, MPI_Wtime();
    double T3, TT3, T4, TT4, T5, TT5, T6, TT6, T7, TT7, T8, TT8, T10, TT10;
    double T0_ave, T1_ave, T2_ave, T3_ave, T4_ave, T5_ave, T10_ave;
    double T6_ave, T7_ave, T8_ave, T9_ave, np;
    double T0a, T1a, T2a, T3a, T4a, T5a, T6a, T7a, T8a, T9a;
    int lszxi, ixi, ix, count_cvg=0; //voro flag: include voronoi tessellation
    double xi_cvg[2], E_elas[2], E_elas_buf[2], E_elas_test;
    int Elascount;
    FILE *of9, *inf1;
    char outfile9[100], output9[100];
    double dxi;
    
    FILE *of2, *of3, *of8, *inf0;
    char outfile[100], output[100], outfile3[100], output3[100], outfile8[100], output8[100];
    
    double **sigma, ***eps;
    double **avesigma, **avepsd, **avepst, **aveps, **ave_eps, **ave_epst, **ave_sigma, **ave_epsd, **avepsts;
    double **xn, **xb, *tau;
    
    fftwnd_mpi_plan plan, iplan;
    fftw_complex *data_fftw, *work, *data_real, *data2_real, *temp_data, *data_core, *data_strain, *data_disp, *work_strain;
    
    //subroutine declarations
    void material(int mflag, double *a, double *mu, double *young, double *c0, double *c1, double *c2, double *c3, double *c4, double *a1, double *a3, double *isf, double *usf, double *nu, double *ll, double *C44, double *C12, double *C11, double *S44, double *S12, double *S11, double C[3][3][3][3]);
    void frec (double *, double *, double *, double, double,double, int, int,int N1, int N2, int N3);
    void setfcc(double **xn, double **xb, int rank);
    void set2D(double **xn, double **xb, int rank, int oflag);
    void set3D1pl(double **xn, double **xb, int rank, int oflag);
    void set3D2sys(double **xn, double **xb, int rank);
    void resolSS( double **sigma, double *tau, double ***eps, int rank, int ND, int NS, double C[][3][3][3], double mu);
    void Bmatrix( double *BB, double *fx, double *fy, double *fz, double **xn, double **xb, double ***eps, double , double , double , double C[][3][3][3], double mu, double nu, double, double, int, int rank, int lsize);
    void Fmatrix(double *FF, double *DD, double *UU, double *fx, double *fy, double *fz, double ***eps, double d1, double d2, double d3, double C[][3][3][3], double mu, double nu, int lN1, int rank, int lsize); /*needs to be called after function Bmatrix*/
    void initial_sxtal(int rank, int lN1, int lxs, double *xi, double *xo, fftw_complex *data_fftw, double *xi_sum, double obsden, double *dslip, double b, int lsize, double nu);
    void rho(double *xi, int lsize, int lN1, int rank);//dislocation density
    void core_ener(int cflag, int it, int rank, int lxs, double An, int lN1, double c0, double c1, double c2, double c3, double c4, double a1, double a3, double **xn, double **xb, fftw_complex *data_fftw, fftw_complex *data_core, double dslip, double b, double *fcore, double *df1core, double *df2core, double *df3core, double *dE_core, double E_core, int itp, double pi, int N1, int N2, int N3, int ND, int NS, int NP, int NT, int NSI, double Bn, double mu, int nsize, int lsize, double size, double isf, char *scratch, double usf);
    void avestrain(double **avepsd, double **avepst, double ***eps, double *xi, int nsize, int lsize, double **sigma, double S11, double S12, double S44, double mu, int lN1, double **ave_epsd, int rank, double **avepsts, int N1, int N2, int N3, int ND, int NS);
    void strain(fftw_complex *data_strain, fftw_complex *data_disp, double *data_eps, fftw_complex *data_fftw, double *FF, double *UU, int tmstep, double **avepst, int lxs, int lN1, int nsize, int lsize, fftw_complex *work_strain, int rank);
    void stress(double *data_epsd, double *data_sigma, fftw_complex *data_strain, double *data_eps, double *xi, double ***eps, double C11, double C12, double C44, int tmstep, double **avesigma, double b, int lxs, int lN1, int nsize, int lsize, int rank, double **ave_sigma, double **sigma, double **avepsts, double *E_elas_test); /*called after function strain*/
    void init_genrand(unsigned long s);
    double genrand_real3(void);
    void verification(int test, double E_elas, double E_elas_test);
    
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    if(rank == 0){
        printf("rank %d, numprocs %d\n", rank, numprocs);
    }
    
    
    if (test == 1){
        inputVariables("PFDD_input_test1.dat");
        if(rank == 0)
            printf("Test stress-strain curve\n----------------\n");
        FILE *infStandard;
        infStandard = fopen(StandardFile,"r");
        if(rank == 0 && infStandard == NULL){
            printf("ERROR IN THE STANDARD FILE, DON'T EXIST OR IS NOT NAMED %s \n", StandardFile);
            exit(1);
        }
        fclose(infStandard);
    }
    
    else if(test == 2){
        if(rank == 0)
            printf("Test elastic energy of screw dipoles\n----------------\n");
        inputVariables("PFDD_input_test2.dat");
    }
    
    T0 = MPI_Wtime();
    
    /*create plan and iplan for fftw*/
    
    plan = fftw3d_mpi_create_plan(MPI_COMM_WORLD, N1, N2, N3, FFTW_FORWARD, FFTW_ESTIMATE);
    
    iplan = fftw3d_mpi_create_plan(MPI_COMM_WORLD, N1, N2, N3, FFTW_BACKWARD, FFTW_ESTIMATE);
    
    /* slab decomposition*/
    
    fftwnd_mpi_local_sizes(plan, &lN1, &lxs, &lnyt, &lyst, &lsize);
    np = numprocs;
    if(rank == 0){
        printf("lN1 %d, lxs %d, lnyt %d, lyst %d, lsize %d\n", lN1, lxs, lnyt, lyst, lsize);
        printf("N1 %d, N2 %d, N3 %d, NT %d, NSI %d, NS %d, ND %d, NP %d, rsm %d, iflag %d, nloops %d, config %d, radius2 %d, mflag %d, cflag %d, An %lf, Bn %lf, oflag %d, setsigma %lf, sigstep %lf, seteps%lf, epstep %lf, voro %d\n",N1, N2, N3, NT, NSI, NS, ND, NP, rsm, iflag, nloops, config, radius2, mflag, cflag, An, Bn, oflag, setsigma, sigstep, seteps, epstep, voro);
    }
    
    /*malloc data vectors*/
    
    data_fftw = (fftw_complex*) malloc(sizeof(fftw_complex) * lsize*NS);
    data_real = (fftw_complex*) malloc(sizeof(fftw_complex) * lsize*NS);
    data2_real = (fftw_complex*) malloc(sizeof(fftw_complex) * lsize*NS);
    temp_data = (fftw_complex*) malloc(sizeof(fftw_complex) * lsize*NS);
    work = (fftw_complex*) malloc(sizeof(fftw_complex) * lsize);
    work_strain = (fftw_complex*) malloc(sizeof(fftw_complex) * lsize);
    data_core = (fftw_complex*) malloc(sizeof(fftw_complex) * lsize*NS);
    data_strain = (fftw_complex*) malloc(sizeof(fftw_complex) * lsize*ND*ND);
    data_disp = (fftw_complex*) malloc(sizeof(fftw_complex)*lsize*ND);
    xi = (double*) malloc(2*(NS)*lsize*sizeof(double));
    xo = (double*) malloc((NS)*lsize*sizeof(double));
    xi_sum = (double*) malloc(2*lsize*sizeof(double));
    fx = (double*) malloc(lsize*sizeof(double));
    fy = (double*) malloc(lsize*sizeof(double));
    fz = (double*) malloc(lsize*sizeof(double));
    f = (double*) malloc((NS)*lsize*sizeof(double));
    r = (double*) malloc((ND)*sizeof(double));
    BB = (double*) malloc((NS)*(NS)*lsize*sizeof(double));
    FF = (double*) malloc((NS)*lsize*(ND)*(ND)*sizeof(double));
    DD = (double*) malloc((NS)*lsize*(ND)*(ND)*sizeof(double));
    UU = (double*) malloc(NS*lsize*ND*sizeof(double));
    fcore = (double*) malloc((NP)*lsize*sizeof(double));
    df1core = (double*) malloc((NP)*lsize*sizeof(double));
    df2core = (double*) malloc((NP)*lsize*sizeof(double));
    df3core = (double*) malloc((NP)*lsize*sizeof(double));
    dE_core = (double*) malloc((NS)*lsize*sizeof(double));
    data_eps = (double*) malloc(2*((ND)*(ND)*lsize)*sizeof(double));
    data_epsd = (double*) malloc(2*((ND)*(ND)*lsize)*sizeof(double));
    data_sigma = (double*) malloc(2*((ND)*(ND)*lsize)*sizeof(double));
    tau = (double*) malloc((NS)*sizeof(double));
    
    xn = (double**) malloc((NS)*sizeof(double));
    xb = (double**) malloc((NS)*sizeof(double));
    
    sigma = (double**) malloc((ND)*sizeof(double));
    epsilon = (double**)malloc((ND)*sizeof(double));
    avesigma = (double**) malloc((ND)*sizeof(double));
    avepsd = (double**) malloc((ND)*sizeof(double));
    avepst = (double**) malloc((ND)*sizeof(double));
    aveps = (double**) malloc((ND)*sizeof(double));
    ave_eps = (double**) malloc((ND)*sizeof(double));
    ave_epst = (double**) malloc((ND)*sizeof(double));
    ave_sigma = (double**) malloc((ND)*sizeof(double));
    ave_epsd = (double**) malloc((ND)*sizeof(double));
    avepsts = (double**) malloc((ND)*sizeof(double));
    
    eps =  (double***) malloc((NS)*sizeof(double));
    
    for(b1=0;b1<NS;b1++){
        xn[b1]=(double*) malloc((ND)*sizeof(double));
        xb[b1]=(double*) malloc((ND)*sizeof(double));
        eps[b1] = (double**) malloc((ND)*sizeof(double));
        for(b2=0;b2<ND;b2++){
            eps[b1][b2] = (double*) malloc((ND)*sizeof(double));
        }
    }
    
    for(b1=0;b1<ND;b1++){
        sigma[b1] = (double*) malloc((ND)*sizeof(double));
        epsilon[b1] = (double*)malloc(ND*sizeof(double));
        avesigma[b1] = (double*) malloc((ND)*sizeof(double));
        avepsd[b1] = (double*) malloc((ND)*sizeof(double));
        avepst[b1] = (double*) malloc((ND)*sizeof(double));
        aveps[b1] = (double*) malloc((ND)*sizeof(double));
        ave_eps[b1] = (double*) malloc((ND)*sizeof(double));
        ave_epst[b1] = (double*) malloc((ND)*sizeof(double));
        ave_sigma[b1] = (double*) malloc((ND)*sizeof(double));
        ave_epsd[b1] = (double*) malloc((ND)*sizeof(double));
        avepsts[b1] = (double*) malloc((ND)*sizeof(double));
    }
    
    //define constants
    L = (double)(N3);
    seed = 243578 + lxs;
    pi = 3.141592654;
    init_genrand(seed);  /* seed the random number generator */
    
    /* Material constants*/
    
    //mflag = 1;
    material(mflag, &a, &mu, &young, &c0, &c1, &c2, &c3, &c4, &a1, &a3, &isf, &usf, &nu, &ll, &C44, &C12, &C11, &S44, &S12, &S11, C);
    b = (a)/sqrt(2.0); //meters-->a is in meters
    
    if(rank == 0){
        printf("mflag %d, a %lf, mu %lf, young %lf, nu %lf, ll %lf, C11 %lf, C12 %lf, C44 %lf, S11 %lf, S12 %lf, S44 %lf, isf %lf, usf %lf\n", mflag, a, mu, young, nu, ll, C11, C12, C44, S11, S12, S44, isf, usf);
    }
    
    nsize = N1*N2*N3;
    lszxi = 2*NS*lsize;
    d1 = 1.0; //size/N1;  //need N1 not lN1 here
    d2 = 1.0; //size/N2;  //in units of b so frequencies are normalized
    d3 = 1.0; //size/N3;
    fflag = 1; //determines whether multiple (fflag=1) or single (fflag!=1 ffts are taken, two different ways to the same result.
    //cflag = 1;
    
    
    if(rank == 0){
        printf("N1 %d, d1 %lf, d2 %lf, d3 %lf, cflag %d\n", N1, d1, d2, d3, cflag);
    }
    
    CD = 1.0; /*dislocation mobility coefficient*/
    dtt = 0.01; /*time increment*/
    
    /*Set the slip systems, frequency and BB matrices */
    
    if(NS == 1 && NP == 1){
        set2D(xn, xb, rank, oflag); //1 slip system edge and screw
    }
    else if(NS == 2 && NP == 2){
        set3D2sys(xn, xb, rank);  //2 slip systems
    }
    else if(NS == 3 && NP == 1){
        set3D1pl(xn, xb, rank, oflag); //3 slip systems on 1 plane edge and screw
    }
    else if(NS == 12 && NP == 4){
        setfcc(xn,xb, rank); //12 slip systems
    }
    else{
        if(rank == 0){
            printf("Direction vectors for the Burger's vectors and slip plane normals have not been included in a subroutine for this number of slip systems or slip planes.\n");
        }
    }
    
    if(rank == 0){
        printf("Calling Frequency Subroutine.\n");
    }
    
    frec(fx, fy, fz, d1, d2, d3, lN1, lxs, N1, N2, N3);
    
    /*Initial Data*/
    
    if(rank == 0){
        printf("Setting Initial Conditions\n");
    }
    
    T3 = MPI_Wtime();
    obsden = 0.1;  //if iflag == 9 this is the dislocation density
    initial_sxtal(rank, lN1, lxs, xi, xo, data_fftw, xi_sum, obsden, &dslip, b, lsize, nu);
    
    //rho(xi, lsize, lN1, rank);
    
    /*fread xitmp: resume the running*/
    if (rsm == 1) {
        static char F[] = "xitmp_P%03d.dat";
        char f[sizeof F +3];
        sprintf(f, F, rank);
        inf1 = fopen(f, "r");
        if (!inf1) {
            printf("File 'xitmp.dat' could not be opened\n");
            perror("xitmp.dat");
            exit(1);
        }
        fread(xi, lszxi, sizeof(double), inf1);
        fclose(inf1);
    }
    else if(rsm == 2){
        char buf[100], buf0[100], buf1[100];
        for(isa=0;isa<NS;isa++){
            sprintf(outfile, "output_it%08d_NS%02d_P%03d.dat", it-1+itp*NT, isa, rank);
            of2 = fopen(outfile, "r");
            if(!of2){
                printf("File %s could not be opened\n", outfile);
                perror("last step output file");
                exit(1);
            }
            if(rank==0){
                fgets(buf, 100, of2);
            }
            for(i=0;i<lsize;i++){
                na0 = 2*i;
                na1 = na0 +1;
                fscanf(of2, "%s %s %s %s %s %s", buf, buf, buf, buf0, buf1, buf);
                xi[na0] = atof(buf0);
                xi[na1] = atof(buf1);
            } //lsize
            
            fclose(of2);
            
        } //isa
    }//rsm==2
    
    MPI_Barrier(MPI_COMM_WORLD);
    
    TT3 = MPI_Wtime() - T3;
    
    
    if(rank == 0){
        printf("Setting interaction matrix.\n");
    }
    
    T8 = MPI_Wtime();
    Bmatrix(BB, fx, fy, fz, xn, xb, eps, d1,  d2 , d3,  C, mu, nu, b, dslip, lN1, rank, lsize); //for calculate dE/dxi
    TT8 = MPI_Wtime() - T8;
    Fmatrix(FF, DD, UU, fx, fy, fz, eps, d1,  d2 , d3,  C, mu, nu, lN1, rank, lsize); //for calculate strain
    
    
    /*Time Step Loop */
    
    if(rank == 0){
        printf("Beginning Time Evolution\n");
    }
    T4 = MPI_Wtime();
    TT1 = 0.0;
    TT2 = 0.0;
    TT5 = 0.0;
    TT6 = 0.0;
    TT7 = 0.0;
    TT9 = 0.0;
    TT10 = 0.0;
    
    for(;itp<NSI;itp++){	//itp may start from any value
        /*Set Applied Stress at stress-controlled condition*/
        
        sigma[0][0]= 0.0;
        sigma[0][1]= 0.0;
        sigma[0][2]= (setsigma + itp*sigstep);
        sigma[1][0]= sigma[0][1];
        sigma[1][1]= 0.0;
        sigma[1][2]= 0.0;
        sigma[2][0]= sigma[0][2];
        sigma[2][1]= sigma[1][2];
        sigma[2][2]= 0.0;
        
        
        if(load == 1){  //strain-controlled case
            epsilon[0][0]= 0.0;
            epsilon[0][1]= 0.0;
            epsilon[0][2]= (seteps + itp*epstep);
            epsilon[1][0]= epsilon[0][1];
            epsilon[1][1]= 0.0;
            epsilon[1][2]= 0.0;
            epsilon[2][0]= epsilon[0][2];
            epsilon[2][1]= epsilon[1][2];
            epsilon[2][2]= 0.0;
        }
        
        xi_cvg[0] = 1.0; xi_cvg[1] = 2.0;
        count_cvg = 0;
        
        if(rank == 0){
            if(load == 0)
                printf("set applied stress sigma[0][2] %lf sigma[2][0] %lf\n", sigma[0][2], sigma[2][0]);
            
            else if(load == 1)
                printf("set strain epsilon[0][2] %lf epsilon[2][0] %lf\n", epsilon[0][2], epsilon[2][0]);
        }
        
        resolSS(sigma, tau, eps, rank, ND, NS, C, mu);
        
        /*Enter 2nd Time Step Loop*/
        
        /* data goes to fftw and is repalced by its fft
         xi is always in real space and is updated every step*/
        
        for(it=0;it<NT;it++){
            tmstep = it + (itp*NT);
            if(rank==0 && it%10 == 0){
                printf("time step, %d %d %d ", it, itp, tmstep);
            }
            
            /* Initialize values each time step*/
            
            for(i=0;i<lsize*2;i++){
                xi_sum[i] = 0.0;
            }
            
            E_core = 0.0;  //Core Energy for each time step.
            E_elas[0] = 0.0; E_elas[1] = 0.0;
            E_elas_buf[0] = 0.0; E_elas_buf[1] = 0.0;
            
            for(i=0;i<lsize*NS;i++){
                dE_core[i] = 0.0;
            }
            
            for(i=0;i<lsize*NP;i++){
                fcore[i] = 0.0;
                df1core[i] = 0.0;
                df2core[i] = 0.0;
                df3core[i] = 0.0;
            }
            
            for(i=0;i<ND;i++)
                for(j=0;j<ND;j++){
                    avesigma[i][j] = 0.0;
                    avepst[i][j] = 0.0;
                    avepsts[i][j] = 0.0;
                    avepsd[i][j] = 0.0;
                    aveps[i][j] = 0.0;
                    ave_eps[i][j] = 0.0;
                    ave_sigma[i][j] = 0.0;
                    ave_epsd[i][j] = 0.0;
                }
            
            //Calculate Average Strain
            
            avestrain(avepsd, avepst, eps, xi, nsize, lsize, sigma, S11, S12, S44, mu, lN1, ave_epsd, rank, avepsts, N1, N2, N3, ND, NS);
            
            //Calculate Core Energy
            T9 = MPI_Wtime();
            
            core_ener(cflag, it, rank, lxs, An, lN1, c0, c1, c2, c3, c4, a1, a3, xn, xb, data_fftw, data_core, dslip, b, fcore, df1core, df2core, df3core, dE_core, E_core, itp, pi, N1, N2, N3, ND, NS, NP, NT, NSI, Bn, mu, nsize, lsize, size, isf, scratch,usf);
            TT9 += MPI_Wtime()-T9;
            //Take forward FFT
            
            if(fflag==1){
                T1 = MPI_Wtime();
                for(isa=0;isa<NS;isa++){
                    psys = isa*lsize;
                    fftwnd_mpi(plan, 1, data_fftw+psys, work, FFTW_NORMAL_ORDER);  /* Forward FFT (multiple)*/
                    //fftwnd_mpi(plan, 1, data_fftw+psys, work, FFTW_TRANSPOSED_ORDER);
                }
                TT1 += (MPI_Wtime() - T1);
            }
            else{
                //T9 = MPI_Wtime();
                for(isa=0;isa<NS;isa++){
                    for(ii=0;ii<lsize;ii++){
                        index = isa + ii*NS;
                        index2 = ii + isa*lsize;
                        temp_data[index] = data_fftw[index2];
                    }
                }
                //TT9 += (MPI_Wtime() - T9);
                
                T1 = MPI_Wtime();
                fftwnd_mpi(plan,NS,temp_data,work,FFTW_NORMAL_ORDER); /*Forward FFT (single)*/
                TT1 += (MPI_Wtime() - T1);
                
                //T9 = MPI_Wtime();
                for(isa=0;isa<NS;isa++){
                    for(ii=0;ii<lsize;ii++){
                        index = isa + ii*NS;
                        index2 = ii + isa*lsize;
                        data_fftw[index2] = temp_data[index];
                    }
                }
                //TT9 += (MPI_Wtime() - T9);
            }
            
            for (i=0; i<lsize*NS; i++){
                data_real[i].re = 0;
                data_real[i].im = 0;
            }
            
            /*Multiply by Interaction Matrix (B-Matrix)*/
            
            T5 = MPI_Wtime();
            for(isa=0;isa<NS;isa++){
                for(isb=0;isb<NS;isb++){
                    for(i=0;i<lN1;i++)
                        for(j=0;j<N2;j++)
                            for(k=0;k<N3;k++){
                                index = k + j*N3 + i*N2*N3;
                                index1  = index + isa*lsize;
                                index2 = index + isb*lsize;
                                nb     = index + isa*lsize + isb*lsize*NS;
                                data_real[index1].re += data_fftw[index2].re * BB[nb];
                                data_real[index1].im += data_fftw[index2].im * BB[nb];
                                E_elas_buf[0] += (data_fftw[index1].re*data_fftw[index2].re + data_fftw[index1].im*data_fftw[index2].im)*BB[nb];
                                E_elas_buf[1] += (data_fftw[index1].im*data_fftw[index2].re - data_fftw[index1].re*data_fftw[index2].im)*BB[nb];
                            }
                }
            }
            
            TT5 += (MPI_Wtime() -T5);
            
            MPI_Allreduce(E_elas_buf, E_elas, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            if(rank ==0){
                E_elas[0] = E_elas[0]*0.5*C44/nsize; //    /N1*b*b;
                for(i=0;i<ND;i++)
                    for(j=0;j<ND;j++){
                        for(k=0;k<ND;k++)
                            for(m=0;m<ND;m++){
                                if(load == 1)
                                    E_elas[0] += nsize*C[i][j][k][m]*(0.5*epsilon[i][j]*epsilon[k][m] - epsilon[i][j]*ave_epsd[k][m]);
                                else if(load == 0){
                                    E_elas[0] += nsize*C[i][j][k][m]*(0.5*avepst[i][j]*avepst[k][m] - avepst[i][j]*ave_epsd[k][m]);
                                } //load=0
                            }//m
                        if(load==0){
                            E_elas[0] -= nsize*mu*sigma[i][j]*avepst[i][j];
                        }
                        
                    }//j
                
                E_elas[0] = E_elas[0]*b*b/N1;
            } //rank==0
            
            /*Strain & Stress Calculations*/
            
            if(it == NT-1){
                if(rank == 0){   //stress-strain curve
                    strcpy(outfile8, scratch);
                    ic = sprintf(output8, "/outsscurve_it%08d_P%03d.dat", tmstep, rank);
                    strcat(outfile8, output8);
                    for(ib=0; ib<strlen(outfile8); ib++){
                        if(outfile8[ib]==' ') outfile8[ib]='0';
                    }
                    
                    of8 = fopen(outfile8,"w");
                }
                
                if(comment == 0){
                    strain(data_strain, data_disp, data_eps, data_fftw, FF, UU, tmstep, avepst, lxs, lN1, nsize, lsize, work_strain, rank);
                    stress(data_epsd, data_sigma, data_strain, data_eps, xi, eps, C11, C12, C44, tmstep, avesigma, b, lxs, lN1, nsize, lsize, rank, ave_sigma, sigma, avepsts, &E_elas_test);
                }
                MPI_Barrier(MPI_COMM_WORLD);
                
                /*average strain*/
                
                for(ida=0;ida<ND;ida++)
                    for (idb=0;idb<ND;idb++){
                        for(i=0;i<lN1;i++)
                            for(j=0;j<N2;j++)
                                for(k=0;k<N3;k++){
                                    na0 = 2*(k + j*N3 + i*N2*N3 + ida*lsize + idb*lsize*ND);
                                    aveps[ida][idb] += data_eps[na0]; //aveps only appears here to get total average strain
                                }
                        
                        MPI_Reduce(&aveps[ida][idb], &ave_eps[ida][idb], 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
                        
                        if(rank == 0){
                            ave_eps[ida][idb] /= nsize;
                            if(load==0){
                                if(ida == 0 && idb == 2){
                                    printf("average strain %lf S*sigma+epsilon0bar %lf\n", ave_eps[ida][idb], avepst[ida][idb]);
                                }
                                fprintf(of8,"%lf %lf %lf %lf ", ave_eps[ida][idb], avepst[ida][idb], ave_sigma[ida][idb]/mu, ave_sigma[ida][idb]);
                            }
                        }
                    }
                
                if(rank==0 && load==0){
                    fprintf(of8, "\n");
                }
                
                else if(load == 1){
                    if(rank == 0){
                        fprintf(of8, "%lf %lf\n", epsilon[0][2], 2*mu*(epsilon[0][2]-ave_epsd[0][2]));
                    }
                }
            }//it==NT-1
            
            //Take inverse FFT
            
            if(fflag==1){
                T2 = MPI_Wtime();
                for(isa=0;isa<NS;isa++){
                    psys = isa*lsize;
                    fftwnd_mpi(iplan, 1, data_real+psys, work, FFTW_NORMAL_ORDER); /* Inverse FFT (multiple)*/
                    //fftwnd_mpi(iplan, 1, data_real+psys, work, FFTW_TRANSPOSED_ORDER);
                }
                TT2 += (MPI_Wtime() - T2);
            }
            else{
                //T9 = MPI_Wtime();
                for(isa=0;isa<NS;isa++){
                    for(ii=0;ii<lsize;ii++){
                        index = isa + ii*NS;
                        index2 = ii + isa*lsize;
                        temp_data[index] = data_real[index2];
                    }
                }
                //TT9 += (MPI_Wtime() - T9);
                
                T2 = MPI_Wtime();
                fftwnd_mpi(iplan,NS,temp_data,work,FFTW_NORMAL_ORDER); /*Inverse FFT (single)*/
                TT2 += (MPI_Wtime() - T2);
                
                //T9 = MPI_Wtime();
                for(isa=0;isa<NS;isa++){
                    for(ii=0;ii<lsize;ii++){
                        index = isa + ii*NS;
                        index2 = ii + isa*lsize;
                        data_real[index2] = temp_data[index];
                    }
                }
                //TT9 += (MPI_Wtime() - T9);
            }
            
            for(isa=0;isa<NS;isa++){
                
                if(it==NT-1||(it%100000==0)){
                    strcpy(outfile,scratch);
                    ic = sprintf(output, "/output_it%08d_NS%02d_P%03d.dat", tmstep, isa, rank);
                    strcat(outfile, output);
                    for(ib=0; ib<strlen(outfile); ib++){
                        if(outfile[ib]==' ') outfile[ib]='0';
                    }
                    of2 = fopen(outfile,"w");
                    if(rank==0){
                      //  fprintf(of2,"zone   I = %d J = %d K = %d\n", N1, N2, N3);
                        fprintf(of2,"zone   I = %d J = %d\n", N1, N2);
                    }
                }
                
                if(rank==0)
                    if(it%dump_3 == 0 && it == NT-1){
                        strcpy(outfile3,scratch);
                        ih = sprintf(output3, "/outaverage_it%08d_NS%02d.dat", tmstep, isa);
                        strcat(outfile3, output3);
                        for(ij=0; ij<strlen(outfile3); ij++){
                            if(outfile3[ij]==' ') outfile3[ij]='0';
                        }
                        
                        of3 = fopen(outfile3,"w");
                    }
                
                xi_ave = 0.0;
                
                T6 = MPI_Wtime();
                
                for(i=0;i<lN1;i++)
                    for(j=0;j<N2;j++)
                        for(k=0;k<N3;k++){
                            index = k + j*N3 + i*N2*N3;
                            index1 = index + isa*lsize;
                            na = 2*index;
                            na0 = 2*index1;
                            na1 = na0+1;
                            
                            //Ginzburg-Landau Equation for real and imag parts
                            
                            if(xo[index1] == 0.0){
                                xi[na0] = xi[na0]-((CD*dtt)*(data_real[index1].re/(nsize) - tau[isa] + dE_core[index1]));
                                xi[na1] = xi[na1]-((CD*dtt)*(data_real[index1].im/(nsize)));
                                xi_sum[na] += xi[na0];
                                xi_sum[na+1] += xi[na1];
                            }
                            
                            data_fftw[index1].re = xi[na0];
                            data_fftw[index1].im = xi[na1];
                            xi_ave += xi[na0];
                            
                            T7 = MPI_Wtime();
                            
                           if(it==NT-1||(it%100000==0)){//give output every itp
 //fprintf(of2, "%d %d %d %e %e %e \n",lxs+i, j,k, xi[na0], xi[na1], xi_sum[na]);
                                if (k==0) {
                                    fprintf(of2,"%d  %d         %e\n",lxs+i,j,xi[na0]);
                                    if (j==N2/2) {
                                        //fprintf(of2,"%d           %e\n",lxs+i,xi[na0]);
                                    }
                                }
                            }
                            TT7 += (MPI_Wtime() - T7);
                        } /*end i,j,k*/
                
                if(it==NT-1||it%100000==0){
                  fclose(of2);
                }
                
                TT6 += (MPI_Wtime() - T6);
                MPI_Reduce(&xi_ave, &xiave, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
                
                if(rank == 0){
                    xiave /= (N1*N2*N3);
                    
                    if(it%10==0 && it<NT-1)
                        printf("average xi %e\n", xiave);
                    
                    if(it%dump_3 == 0 && it == NT-1){
                        if(isa<NS-1)
                            fprintf(of3, "%d %lf\n", tmstep, xiave);
                        else
                            fprintf(of3, "%d %lf %e\n", tmstep, xiave, E_elas[0]);
                        
                        fclose(of3);
                    }
                    
                    if(isa == NS/3){
                        xi_cvg[0] = xi_cvg[1];
                        xi_cvg[1] = xiave;
                    } //isa==NS/3
                    xiave = 0.0;
                } //rank==0
                
                MPI_Bcast(xi_cvg, 2, MPI_DOUBLE, 0, MPI_COMM_WORLD);
                
                
            } /* end isa */
            
            if(rank == 0 && it == NT-1){
                fclose(of8);
                printf("E_elas = %e Im %e\n", E_elas[0], E_elas[1]);
            }
            
            if(itp%80 == 0 && (it == NT-1) && itp>0){		//binery check pointing
                strcpy(outfile9,scratch);
                ixi = sprintf(output9, "/xitmp_it%08d_P%03d.dat", tmstep, rank);
                strcat(outfile9, output9);
                for(ix=0; ix<strlen(outfile9); ix++){
                    if(outfile9[ix]==' ') outfile9[ix]='0';
                }
                of9 = fopen(outfile9,"w");
                fwrite(xi, lszxi, sizeof(double), of9);
                fclose(of9);
            }
            if(it < NT-1)
                if(iflag != 12){
                    if(fabs((xi_cvg[1]-xi_cvg[0])/xi_cvg[1]) < 1.0E-20){   // check the convergence of phase field xi
                        count_cvg += 1;
                        if(count_cvg > 10){
                            count_cvg = 0;
                            if(rank == 0)
                                printf("phase field converged at time step %d \n", it+itp*NT);
                            it = NT -2;
                        }
                    }
                    else
                        count_cvg = 0;   // reset count_cvg
                }
                else{
                    it = NT -2;
                }
            
            TT6 += (MPI_Wtime() - T6);
        }  /*end it*/
        
        /*calculate dislocation density*/
        //rho(xi, lsize, lN1, rank);
        
    } /*end itp*/
    TT4 = MPI_Wtime() - T4;
    /*free malloc'ed memory*/
    if(rank == 0)
        if(test > 0){
            if(comment > 0){
                printf("stress subroutine is not called, E_elas_test is not calculated \n");
                exit(1);
            }
            verification(test, E_elas[0], E_elas_test);
        }
    
    free(data_fftw);
    free(data_eps);
    free(data_epsd);
    free(data_sigma);
    free(work);
    free(work_strain);
    free(data_real);
    free(data2_real);
    free(temp_data);
    free(xi);
    free(xo);
    free(xi_sum);
    free(fx);
    free(fy);
    free(fz);
    free(f);
    free(r);
    free(BB);
    free(FF);
    free(DD);
    free(UU);
    free(data_core);
    free(data_strain);
    free(data_disp);
    free(fcore);
    free(df1core);
    free(df2core);
    free(df3core);
    free(dE_core);
    free(tau);
    
    for(b1=0;b1<NS;b1++){
        free(xn[b1]);
        free(xb[b1]);
        for(b2=0;b2<ND;b2++){
            free(eps[b1][b2]);
        }
    }
    
    for(b1=0; b1<NS; b1++) {
        free(eps[b1]);
    }
    
    for(b1=0;b1<ND;b1++){
        free(sigma[b1]);
        free(avesigma[b1]);
        free(ave_sigma[b1]);
        free(avepsd[b1]);
        free(avepst[b1]);
        free(ave_epsd[b1]);
        free(avepsts[b1]);
        free(aveps[b1]);
        free(ave_eps[b1]);
        free(ave_epst[b1]);
    }
    
    free(xn);
    free(xb);
    
    free(sigma);
    free(avesigma);
    free(ave_sigma);
    free(avepsd);
    free(avepst);
    free(ave_epsd);
    free(avepsts);
    free(aveps);
    free(ave_eps);
    free(ave_epst);
    
    free(eps);
    
    fftwnd_mpi_destroy_plan(plan);
    fftwnd_mpi_destroy_plan(iplan);
    
    TT0 = MPI_Wtime() - T0;
    
    MPI_Barrier(MPI_COMM_WORLD);
    
    MPI_Reduce(&TT0, &T0_ave, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&TT1, &T1_ave, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&TT2, &T2_ave, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&TT3, &T3_ave, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&TT4, &T4_ave, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&TT5, &T5_ave, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&TT6, &T6_ave, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&TT7, &T7_ave, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&TT8, &T8_ave, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&TT9, &T9_ave, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    
    MPI_Barrier(MPI_COMM_WORLD);
    
    if(rank == 0){
        T0_ave /= np;
        T1_ave /= np;
        T2_ave /= np;
        T3_ave /= np;
        T4_ave /= np;
        T5_ave /= np;
        T6_ave /= np;
        T7_ave /= np;
        T8_ave /= np;
        T9_ave /= np;
        
        printf("Dimension = %d x3 and NP = %f \n", N1, np);
        printf("Time Summary - Average Times\n");
        printf("Time to create B-matrix = %8.3f seconds\n", T8_ave);
        printf("Time to create initial data = %8.3f seconds\n", T3_ave);
        printf("Time for time stepping loop = %8.3f seconds\n", T4_ave);
        printf("Time to take forward FFT = %8.3f seconds\n", T1_ave);
        printf("Time to multiply data by B-matrix = %8.3f seconds\n", T5_ave);
        printf("Time to take inverse FFT = %8.3f seconds\n", T2_ave);
        printf("Time to core_ener = %8.3f seconds\n", T9_ave);
        printf("Time to update data = %8.3f seconds\n", T6_ave);
        printf("Time to write output file = %8.3f seconds\n", T7_ave);
        printf("Time for entire code to run = %8.3f seconds\n", T0_ave);
    }
    
    MPI_Finalize();
    return 0;
}

/********************************************************************/

/* subroutines*/

/*********************************************************************/

void initial_sxtal(int rank, int lN1, int lxs, double *xi, double *xo, fftw_complex *data_fftw, double *xi_sum, double obsden, double *dslip, double b, int lsize, double nu)
{
    
    /*iflag == 1 -> 1 dislocation loop on 1 slip system
     iflag == 2 -> 2 dislocations on 2 slip systems
     iflag == 3 -> 4 obstacles (set for 128x128x128)
     iflag == 4 -> interface (passivated film), random distribution of obstacles
     iflag == 5 -> interface with columnar obstacle lines (like iflag==4)
     iflag == 6 -> interface (passivated film), random dislocation distribution (like iflag==4)
     iflag == 7 -> infinitely long straight unit dislocation on a 111 plane (1 ss) dislocation line parallel to y-axis -> cflag == 1 (stress/strain verification)
     iflag == 8 -> infinitely long straight unit dislocation on a 111 plane (3 ss) dislocation line parallel to y-axis -> cflag == 2 (paritals-verification)
     iflag == 9 -> single grain (3 grain boundaries), (stress/strain curves)
     iflag == 10 -> interface (passivated film) with columnar grain boundaries generated by matlab code -- grid.m
     iflag == 11 -> two screw dislocation (displacement PFDD-MD)
     iflag == 12 -> pure elastic phase field doesn't evolve
     iflag == 13 -> bimodal grains (large grains + small grains)
     iflag == 14 -> 4 planes and 12 slip system active
     iflag == 15 -> Frank-Read source
     iflag == 16 -> 1 dislocation line on 1 slip system w BC
     iflag == 17 -> 3 slip system 1 slip plane, 0-1-0 step, full dislocation splitting to partials
     iflag == 18 -> 1 slip system 1 slip plane, 0-1-0 step, 2 full dislocation annihilate. test elastic part*/
    
    
    void init_genrand(unsigned long s);
    double genrand_real3(void);
    
    int is, i, j, k, im, jm, km, ism, ir, ia, ib, indexgb, grainb, num, count, *gb3d;
    int na0, na, na1, index, index2d, index3d, indexm, *nodes, rtn_val0, rtn_val1, n_obs[7];
    double nlx, nly, c0, c1, a, alpha, zeta, eta, d, *gb2d;
    unsigned long seed;
    FILE *of0, *fgrid, *inf1;
    char infile[100], input[100], c[10];
    for(i=0;i<7;i++) {
        n_obs[i] = 0;
    }
    gb2d = (double*)malloc(N1*N2*sizeof(double));
    gb3d = (int*)malloc(lsize*sizeof(int));
    seed = 243578 + lxs;
    num = 0.0;
    count = 0.0;
    a = 1.0; //approximate core region in iflag == 7
    d = 2.0;
    zeta = (d/(2.0*(1-nu))); //approximate core region in iflag == 7 PN-model edge
    eta = d/2.0;  //approximate core region in iflag == 7 PN-model screw
    alpha = 1.0;
    init_genrand(seed);  /* seed the random number generator */
    
    if(rank == 0){
        printf("Initial data ------------------");
    }
    
    if(iflag == 10){
        nodes = (int*) malloc(N1*N2*sizeof(int));
        if(rank == 0){
            fgrid = fopen("inputnodes.dat", "r");
            if(fgrid == NULL){
                printf("File 'inputnodes.dat' could not be opened\n");
                perror("inputnodes.dat");
                exit(1);
            }
            else{
                while(fgets(c, 10, fgrid) != NULL){
                    //printf("node is: %d \n",atoi(c));
                    while(num != atoi(c) && num <= N1*N2){
                        num++;
                    }
                    if(num == atoi(c)){
                        nodes[count] = num;
                        //printf("node is: %d, %d %d %d\n", atoi(c), num, count, nodes[count]);
                        count++;
                    }
                }
            }
        }
        rtn_val1 = MPI_Bcast(&count, 1, MPI_INT, 0, MPI_COMM_WORLD);
        rtn_val0 = MPI_Bcast(nodes, count, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Barrier(MPI_COMM_WORLD);
        printf("rank %d, rtn_val %d, %d %d %d %d\n", rank, rtn_val0, rtn_val1, count, nodes[0], nodes[count]);
        
        indexgb = 0;
        while(indexgb <= count){
            grainb = nodes[indexgb];
            printf("node: %d, %d, %d, %d, rank %d\n", grainb, nodes[indexgb], indexgb, count, rank);
            //if(num == atoi(c)){
            jm = (grainb-1)/N2;
            im = grainb-jm*N1-1;
            //printf("%d, im=%d, jm=%d \n", grainb, im, jm);
            if(im < N1 && jm < N2 && im >= lxs+i && im < (lxs+i)+lN1){
                for(ism=0;ism<NS;ism++)
                    for(km=0;km<N3;km++)
                    {
                        indexm = (im-lxs)*N2*N3+jm*N3+km+ism*lsize;
                        xo[indexm] = 1.0;
                    }
                //printf("indexm %d, rank %d, im %d, jm %d lsx %d\n", indexm, rank, im, jm, lxs);
            }
            //}
            indexgb++;
        }
    }
    
    if(voro == 1) {
        sprintf(input, "voronoi_%03d.dat", rank);
        inf1 = fopen(input, "r");
        for (i=0; i<lsize; i++) {
            fscanf(inf1, "%d %d %d %d", &is, &is, &is, &gb3d[i]);
        }
        fclose(inf1);
        
        /*if(rank == 0) { //2-D voronoi tessilation
         inf1 = fopen("voronoi.dat", "r");
         if(inf1 == NULL) {
         printf("File 'voronoi.dat' could not be opened\n");
         perror("voronoi.dat");
         exit(1);
         }
         for(i=0;i<N1*N2;i++) {
         fscanf(inf1, "%lf", &gb2d[i]);
         }
         fclose(inf1);
         }*/
    }
    //MPI_Bcast(gb2d, N1*N2, MPI_DOUBLE, 4, MPI_COMM_WORLD);
    
    for(is=0;is<NS;is++){
        /*open input file*/
        strcpy(infile,scratch);
        ia = sprintf(input, "/input_NS%02.0d_P%03.0d.dat", is, rank);
        strcat(infile,input);
        for(ib=0; ib<strlen(infile); ib++){
            if(infile[ib]==' ') infile[ib]='0';
        }
        
        of0 = fopen(infile,"w");
        
        if(rank==0){
            fprintf(of0,"zone  I = %d J = %d K = %d \n",N1, N2, N3);
        }
        
        //printf("is %d\n",is);
        for(i=0;i<lN1;i++)
            for(j=0;j<N2;j++)
                for(k=0;k<N3;k++){
                    na0 = 2*(i*N2*N3+j*N3+k+is*lsize);
                    na = 2*(i*N2*N3+j*N3+k);
                    index = i*N2*N3+j*N3+k+is*lsize;
                    index2d = j + (lxs + i)*N2;
                    index3d = k + j*N3 + i*N2*N3;
                    na1 = na0+1;
                    xi[na0] = 0.0;      /*(k+1)+(j+1)*10.0+(i+1)*100.0*/
                    xi[na1] = 0.0;
                    xo[index] = 0.0;
                    
                    if(k==N3/2){
                        if(iflag == 1){  //1 slip system
                            ir = (lxs+i-N1/2)*(lxs+i-N1/2)+(j-N2/2)*(j-N2/2);
                            if(ir<=N1/2){
                                xi[na0]=1.0;
                                xi_sum[na] = xi_sum[na] + xi[na0];
                            }
                            //write input file
                            fprintf(of0, "%d %d %lf %lf %lf\n", lxs+i, j, xi[na0], xi[na1], xi_sum[na]);
                        } //iflag == 1
                        
                        else if(iflag == 2){  //2 slip system 2 dislocations
                            if(is==0){
                                ir = (lxs+i-N1/4)*(lxs+i-N1/4) + (j-N2/2)*(j-N2/2);
                            }
                            else{
                                ir = (lxs+i-(3*N1/4))*(lxs+i-(3*N1/4)) + (j-N2/2)*(j-N2/2);
                            }
                            if(ir<=N1/4){
                                xi[na0]=1.0;
                                xi_sum[na] += xi[na0];
                            }
                        } //iflag == 2
                        
                        else if(iflag == 3){ //Four obstacles (32x32x32)
                            ir = (lxs+i-N1/2)*(lxs+i-N1/2)+(j-N2/2)*(j-N2/2);
                            if(ir<=N1/2){
                                xi[na0]=1.0;
                                //xi_sum[na] = xi_sum[na] + xi[na0];
                            }
                            if(lxs+i== 16 /*64*/){
                                if(j== 11/*55*/ || j== 21 /*73*/){
                                    xo[index] = 1.0;
                                    xi[na0] = 0.0;
                                }
                            }
                            else if(j== 16 /*64*/){
                                if(lxs+i== 11 /*56*/ || lxs+i== 21 /*72*/){
                                    xi[na0] = 0.0;
                                    xo[index] = 1.0;
                                }
                            }
                        } //iflag == 3
                    } // k == N3/2
                    
                    if(iflag == 4){ //interface with obstacles
                        xo[index] = genrand_real3();
                        //printf("%d %d %d %d %lf\n",i, j, k, index, xo[index]);
                        if(xo[index] <= obsden){
                            xo[index] = 1.0;
                        }
                        else{
                            xo[index] = 0.0;
                        }
                        if(k <= (N3-1)/6 || k >= 5*(N3-1)/6){
                            xo[index] = 1.0;
                        }
                        /*if(j <= (N2-1)/6 || j >= 5*(N2-1)/6){
                         xo[index] = 1.0;
                         }*/
                        
                    } //iflag == 4
                    
                    if(iflag == 5){ //interface with columnar obstacles
                        nlx = (2.0*N1)/((obsden*N1)+1);
                        nlx = floor(nlx + 0.5);
                        nly = (2.0*N2)/((obsden*N2)+1);
                        nly = floor(nly +0.5);
                        //printf("%d %d %d %d %lf\n",i, j, k, index, xo[index]);
                        if(fmod(lxs+i,nlx) == 0 || fmod(j,nly) == 0){
                            xo[index] = 1.0;
                        }
                        else if(lxs+i == (N1-1) || j == (N2-1)){
                            xo[index] = 1.0;
                        }
                        else{
                            xo[index] = 0.0;
                        }
                        if(k <= (N3-1)/6 || k >= 5*(N3-1)/6){
                            xo[index] = 1.0;
                        }
                    } //iflag == 5
                    
                    if(iflag == 6){  //interface with random dislocation distribution
                        if(fmod((double)(k),10.0) == 0.0){
                            xi[na0] = genrand_real3();
                            //printf("%d %d %d %d %lf\n",i, j, k, index, xo[index]);
                            if(xi[na0] <= obsden){
                                xi[na0] = 1.0;
                            }
                            else{
                                xi[na0] = 0.0;
                            }
                        }
                        else{
                            xi[na0] = 0.0;
                        }
                        if(k <= (N2-1)/6 || k >= 5*(N2-1)/6){
                            xo[index] = 1.0;
                            xi[na0] = 0.0;
                        }
                        else{
                            xo[index] = 0.0;
                        }
                    } //iflag == 6
                    
                    if(iflag == 7){ //infinitely long dislocation 1 slip system dislocation line parallel to y-axis, core region approimately by a = 1. -> cflag == 1 (stress/strain-verification)
                        if(k == N3/2){
                            //xi[na0] = 1.0-(1.0/(1+exp(-(a*((lxs+i)-(N1/2))))));
                            //xi[na0] = 0.5+(-(1/pi)*atan(((lxs+i)-(N1/2))/zeta)); //PN-m	\odel edge
                            //xi[na0] = 0.5+(-(1/pi)*atan(((lxs+i)-(N1/2))/eta)); //PN-mo\del screw
                            if((lxs+i) <= N1/2){
                                xi[na0] = 1.0; //step function
                            }
                            
                        }
                    } //iflag == 7
                    
                    if(iflag == 8){  //infinitely long dislocation 3 slip systems dislocation line parallel to y-axis -> cflag == 2 (partials-verification)
                        if(is == 1){
                            //if(j >= N2/2 && k == N3/2){
                            if((lxs+i) <= N1/2 && k == N3/2){
                                //if((lxs+i) > (N1/4)+(1.0/alpha) && (lxs+i) < (3*N1/4)-(1.0/alpha) && j == N2/2){
                                xi[na0] = 1.0;
                                if((lxs+i) == N1/2 && j == 0){
                                    xo[index] = 1.0;
                                }
                                else if((lxs+i) == N1/2 && j == N2-1){
                                    xo[index] = 1.0;
                                }
                            }
                            //else if((lxs+i) >= N1/4 && (lxs+i) <= (N1/4)+(1.0/alpha) && j == N2/2){
                            //xi[na0] = (alpha)*((lxs+i)-(N1/4));
                            //}
                            //else if((lxs+i) >= (3*N1/4)-(1.0/alpha) && (lxs+i) <= (3*N1/4) && j == N2/2){
                            //xi[na0] = 1.0 - (alpha)*((lxs+i)-((3*N1/4)-(1.0/alpha)));
                            //}
                            //else{
                            //xi[na0] = 0.0;
                            //}
                        }
                        else if(is == 0){
                            //if(j >= N2/2 && k == N3/2){
                            //if((lxs+i) >= N1/2 && k == N3/2){
                            xi[na0] = 0.0;
                            /* if((lxs+i) == N1/2 && k == N3/2 && j == 0){
                             xo[index] = 1.0;
                             }
                             else if((lxs+i) == N1/2 && k == N3/2 && j == N2-1){
                             xo[index] = 1.0;
                             }*/
                            //}
                        }
                        else if(is == 2){
                            //if(j >= N2/2 && k == N3/2){
                            //if((lxs+i) >= N1/2 && k == N3/2){
                            xi[na0] = 0.0;
                            //}
                            /*if((lxs+i) == N1/2 && k == N3/2 && j == 0){
                             xo[index] = 1.0;
                             }
                             else if((lxs+i) == N1/2 && k == N3/2 && j == N2-1){
                             xo[index] = 1.0;
                             }*/
                        }
                    } //iflag == 8
                    
                    if(iflag == 9){ //one grain (inpenatrable grain boundary), every 4th plane is active, number of initial loops is variable (stress/strain curves)
                        *dslip = 4.0;
                        if(k%4==0){
                            if(gb3d[index3d] == 1){
                                /**********random obstacles**********/
                                xo[index] = genrand_real3();
                                if(xo[index] <= 0.032){			//random obs density 0.1
                                    n_obs[6]++;
                                    xo[index] = 1.0;
                                    xi[na0] = 1 - rand() % 3;	// -1~1
                                    //xi[na0] = -3 + rand() % 7;		// -3~3
                                }
                            }
                            
                            /**********random obstacles**********/
                            
                            /* count # of obstacles*/
                            /*if(xi[na0] == -3.0) n_obs[0]++;
                             if(xi[na0] == -2.0) n_obs[1]++;
                             if(xi[na0] == -1.0) n_obs[2]++;
                             if(xi[na0] == 1.0) n_obs[3]++;
                             if(xi[na0] == 2.0) n_obs[4]++;
                             if(xi[na0] == 3.0) n_obs[5]++;*/
                            
                            if(nloops == 1){
                                ir = (lxs+i-N1/2)*(lxs+i-N1/2)+(j-N2/2)*(j-N2/2);
                                if(ir<=radius2){
                                    xi[na0] = 1.0;
                                }
                                else{
                                    xi[na0] = 0.0;
                                }
                            }
                            if(nloops >1){
                                for(c0=0;c0<(nloops/2);c0++)
                                    for(c1=0;c1<(nloops/2);c1++){
                                        ir = (lxs+i-((2*c0+1)*N1)/nloops)*(lxs+i-((2*c0+1)*N1)/nloops)+(j-((2*c1+1)*N2)/nloops)*(j-((2*c1+1)*N2)/nloops);
                                        if(ir<=radius2){
                                            xi[na0] = 1.0;
                                        }
                                    }
                            }
                            
                            if(gb2d[index2d] == 1.0) {
                                xo[index] = 1.0;
                            }
                        }
                        else{
                            xi[na0] = 0.0;
                            xo[index] = 1.0;
                        }
                        
                    } //iflag == 9
                    
                    if(iflag == 10){ //passivated films with columnar grains--matlab generated obstacle array
                        //if(fmod((double)(k),10.0) == 0.0){
                        xi[na0] = genrand_real3();
                        //printf("%d %d %d %d %lf\n",i, j, k, index, xo[index]);
                        if(xi[na0] <= obsden){
                            xi[na0] = 1.0;
                        }
                        else{
                            xi[na0] = 0.0;
                        }
                        //}
                        if(k <= (N3-1)/6 || k >= 5*(N3-1)/6){
                            xo[index] = 1.0;
                        }
                    } //iflag == 10
                    
                    if(iflag == 11){
                        *dslip = N3;
                        int disl_width=0, nb;	//disl core spread
                        if(k == N3/2){
                            if(is == NS/3){                                  //separated screw
                                if(j>=(N2/4+disl_width/2) && j<(3*N2/4-disl_width/2)){
                                    xi[na0] = 1.0;
                                }
                                for(nb=1;nb<disl_width;nb++){
                                    if(j==(N2/4-disl_width/2 +nb)){
                                        xi[na0] = 1.0*nb/disl_width;
                                    }
                                    if(j==(N2/4*3-disl_width/2-1 +nb)){
                                        xi[na0] = 1-1.0*nb/disl_width;
                                    }
                                }
                            }
                        }
                        
                    }//iflag == 11
                    
                    if(iflag == 12){
                        *dslip = 1.0;
                        xi[na0] = 0.0;
                        xo[index] = 1.0;
                    } // iflag==12 pure elastic
                    
                    if(iflag == 13){  /*********bimodal configurations**********/
                        *dslip = 4.0;
                        if((k-2)%(int)(*dslip)==0){
                            xi[na0] = 0.0;
                            xo[index] = 0.0;
                        }
                        else{
                            xi[na0] = 0.0;
                            xo[index] = 1.0;
                        }
                        if(config == 0) {
                            if(lxs+i == 0 || lxs+i == (N1-1) || j == 0 || j == (N2-1) || k == 0 || k == (N3-1)){
                                //if(lxs+i == 0 || j == 0 || k == 0){
                                xi[na0] = 0.0;
                                xo[index] = 1.0;
                            }
                        }
                        
                        if(config == 1) {
                            if(lxs+i == 0 || lxs+i == (N1-1) || j == 0 || j == (N2-1) || k == 0 || k == (N3-1) || lxs+i == N1/2 || (j == N1/2 && (lxs+i)>N1/2)){
                                xi[na0] = 0.0;
                                xo[index] = 1.0;
                            }
                        }
                        
                        if(config == 2) {
                            if(lxs+i == 0 || lxs+i == (N1-1) || j == 0 || j == (N2-1) || k == 0 || k == (N3-1) || lxs+i == N1/2 || ((j == N1/2 || j == N1/4 || j == N1*3/4) && (lxs+i)>N1/2)){
                                xi[na0] = 0.0;
                                xo[index] = 1.0;
                            }
                        }
                        
                        if(config == 3) {
                            if(lxs+i == 0 || lxs+i == (N1-1) || j == 0 || j == (N2-1) || k == 0 || k == (N3-1) || lxs+i == N1/2 || j == N1/2 || ((j == N1/4 || j == N1*3/4) && (lxs+i)>N1/2)){
                                xi[na0] = 0.0;
                                xo[index] = 1.0;
                            }
                        }
                        
                        if(config == 4) {
                            if(lxs+i == 0 || lxs+i == (N1-1) || j == 0 || j == (N2-1) || k == 0 || k == (N3-1) || lxs+i == N1/2 || lxs+i == N1*3/4 || ((j == N1/2 || j == N1/4 || j == N1*3/4) && (lxs+i)>N1/2)){
                                xi[na0] = 0.0;
                                xo[index] = 1.0;
                            }
                        }
                        
                        if (config == 5) {
                            if (lxs+i == 0 || lxs+i == (N1-1) || j == 0 || j == (N2-1) || k == 0 || k == (N3-1) || lxs+i == N1*3/4 || ((j == N1/2 || j == N1/4 || j == N1*3/4) && lxs+i>N1*3/4)) {
                                xi[na0] = 0.0;
                                xo[index] = 1.0;
                            }
                        }
                        
                    }//iflag==13 bimodal grains
                    
                    if(iflag == 14){
                        if(lxs+i+j+k == N1-1 || lxs+i+j-k == N1-1 || lxs+i-j+k == N1-1 || j+k-lxs-i == N1-1)
                            xo[index] = 0.0;
                        else
                            xo[index] = 1.0;
                    } // iflag == 14
                    
                    if(iflag == 15){
                        *dslip = 1.0;
                        if(k==N3/2){
                            xo[index] = 0.0;
                            if(is==NS/3){
                                if(lxs+i>N1/3 && lxs+i<N1/3*2 && j>N2/3 && j<N2/3*2){
                                    xi[na0] = 1.0;
                                    xo[index] = 1.0;
                                }
                            }
                        }
                        else
                            xo[index] = 1.0;
                    }// iflag == 15
                    
                    if(iflag == 16){
                        *dslip = 1.0;
                        if(k==N3/2){
                            xo[index] = 0.0;
                            if(lxs+i<N1/2){
                                xi[na0] = 1.0;
                            }
                            
                            if(j==0 || j==N2-1){
                                xo[index] = 1.0;
                            }
                        }
                    }//iflag == 16
                    
                    if(iflag == 17){
                        *dslip = 1.0;
                        if (k==N3/2) {
                            xo[index] = 0.0;
                            
                            if (is==1) {
                                if (lxs+i<=N1*3/4&&lxs+i>=N1/4) {
                                    xi[na0] = 1.0;
                                }
                                if (lxs+i==N1/2||lxs+i==0||lxs+i==N1-1) {
                                    xo[index] = 1.0;
                                    
                                }
                                
                            }else if(is==0){
                                xi[na0] = 0.0;
                                if(lxs+i==N1/2||lxs+i==0||lxs+i==N1-1){
                                    xo[index] = 1.0;
                                }
                            }else if(is==2){
                                xi[na0] = 0.0;
                                if(lxs+i==N1/2||lxs+i==0||lxs+i==N1-1){
                                    xo[index] = 1.0;
                                }
                            }
                        }else{
                            xo[index] = 1.0;
                        }
                    }// iflag == 17
                    
                    
                    if(iflag == 18){
                        *dslip = 1.0;
                        if (k==N3/2) {
                            xo[index] = 0.0;
                            
                            if (is==0) {
                                if (lxs+i<=N1*3/4&&lxs+i>=N1/4) {
                                    xi[na0] = 1.0;
                                }
                                if (lxs+i==N1/2) {
                                    xo[index] = 1.0;
                                }
                            }
                        }else{
                            xo[index] = 1.0;
                        }
                    }// iflag == 18
                    
                    
                    
                    if(iflag == 20){
                        *dslip = 1.0;
                        if (k==0) { //on xy plane, k = 0
                            
                            //xo[index] = 0.0;

                          // if (lxs+i <=2/3.*N1 && lxs+i >=1/3.*N1 &&(j >=N2*3/4.|| j <=N2*1/4.)){
                         
                           //if (((lxs+i)-N1*0.45)*((lxs+i)-N1*0.45)+(j-N2*0.)*(j-N2*0.) <= (N2*0.275)*(N2*0.275) || ((lxs+i)-N1*0.45)*((lxs+i)-N1*0.45)+(j-N2*1)*(j-N2*1) <= (N2*0.275)*(N2*0.275)) { //d=0.45N2

                        // if (((lxs+i)-N1*0.45)*((lxs+i)-N1*0.45)+(j-N2*0.)*(j-N2*0.) <= (N2*0.3)*(N2*0.3) || ((lxs+i)-N1*0.45)*((lxs+i)-N1*0.45)+(j-N2*1)*(j-N2*1) <= (N2*0.3)*(N2*0.3)) { //d=0.4N2

                         // if (((lxs+i)-N1*0.45)*((lxs+i)-N1*0.45)+(j-N2*0.)*(j-N2*0.) <= (N2*0.25)*(N2*0.25) || ((lxs+i)-N1*0.45)*((lxs+i)-N1*0.45)+(j-N2*1.)*(j-N2*1.) <= (N2*0.25)*(N2*0.25)) { //d=0.5N2

                        // if (((lxs+i)-N1*0.45)*((lxs+i)-N1*0.45)+(j-N2*0.)*(j-N2*0.) <= (N2*0.325)*(N2*0.325) || ((lxs+i)-N1*0.45)*((lxs+i)-N1*0.45)+(j-N2*1)*(j-N2*1) <= (N2*0.325)*(N2*0.325)) { //d=0.35N2

                        // if (((lxs+i)-N1*0.45)*((lxs+i)-N1*0.45)+(j-N2*0.)*(j-N2*0.) <= (N2*0.35)*(N2*0.35) || ((lxs+i)-N1*0.45)*((lxs+i)-N1*0.45)+(j-N2*1)*(j-N2*1) <= (N2*0.35)*(N2*0.35) || (lxs+i)>(N1-30)) { //d=0.3N2

                           if ((lxs+i)>(N1-30)){
                               xo[index]=1.0;
                            }else{
                               xo[index] = 0.0;
                            }

                            if (is==0) { //slip plane no.1

                                if (lxs+i<=N1*1/4.) {
                                    xi[na0] = 1.0;
                                }
                                if (lxs+i==N1/8.) {
                                    xo[index] = 0;
                                }
                            }
                        }else{
                            xo[index] = 1.0;
                        }
                        
                    }// iflag == 20
                    
                    if(voro == 1){
                        if(gb3d[index3d] == 1)
                            xo[index] = 1.0;
                    }
                    fprintf(of0, "%d %d %d %lf %lf %lf\n", lxs+i, j, k, xi[na0], xi[na1], xo[index]);
                    data_fftw[index].re = xi[na0];
                    data_fftw[index].im = xi[na1];
                    
                } //ijk
    } //is
    
    /*# of fixed obs for xi=-3 -2 -1 1 2 3 and in total*/
    //printf("rank%d %d %d %d %d %d %d %d\n", rank, n_obs[0], n_obs[1], n_obs[2], n_obs[3], n_obs[4], n_obs[5], n_obs[6]);
    MPI_Barrier(MPI_COMM_WORLD);
    
    fclose(of0);
    if(iflag == 10){
        free(nodes);
        if(rank == 0){
            fclose(fgrid);
        }
    }
    
    if(rank == 0){
        printf("Leaving Initial\n");
    }
    
    free(gb2d);
    free(gb3d);
    MPI_Barrier(MPI_COMM_WORLD);
    
    return;
}
/******************************************************************/

void rho(double *xi, int lsize, int lN1, int rank){
    int xleft, xright, yleft, yright, zleft, zright, dx, naxright, naxleft, nayright, nayleft, nazright, nazleft;
    int i, j, k, is, bca;
    double burger, x_r[N1/2], x_l[N1/2], y_r[N1/2], y_l[N1/2], rho[NS], rho_tot;
    burger = 2.49E-10;
    dx = 2;
    rho_tot = 0.0;
    for (i=0; i<NS; i++) {
        rho[i] = 0.0;
    }
    for (i=0; i<N1/2; i++) {
        x_r[i] = 0.0;
        x_l[i] = 0.0;
        y_r[i] = 0.0;
        y_l[i] = 0.0;
    }
    for (is=0; is<NS; is++) {
        for(i=0;i<lN1;i++)
            for(j=0;j<N2;j++)
                for(k=0;k<N3;k++){
                    
                    yright = j+dx/2;
                    yleft = j-dx/2;
                    //zright = k+dx/2;
                    //zleft = k-dx/2;
                    //if(xleft<0){xleft += N1;}
                    if(yleft<0){yleft += N2;}
                    //if(zleft<0){zleft += N3;}
                    //if(xright>=N1){xright += -N1;}
                    if(yright>=N2){yright += -N2;}
                    //if(zright>=N2){zright += -N3;}
                    nayright = 2*(k+yright*N3+i*N2*N3+is*lsize);
                    nayleft = 2*(k+yleft*N3+i*N2*N3+is*lsize);
                    y_r[rank] = xi[nayright];
                    y_l[rank] = xi[nayleft];
                    for (bca=0; bca<N1/2; bca++){
                        MPI_Bcast(&y_r[bca], 1, MPI_DOUBLE, bca, MPI_COMM_WORLD);
                        MPI_Bcast(&y_l[bca], 1, MPI_DOUBLE, bca, MPI_COMM_WORLD);
                    }
                    MPI_Barrier(MPI_COMM_WORLD);
                    //nazright = 2*(zright+j*N3+i*N2*N3+is*N1*N2*N3);
                    //nazleft = 2*(zleft+j*N3+i*N2*N3+is*N1*N2*N3);
                    
                    if (i==0) {
                        xright = i+dx/2;
                        naxright = 2*(k+j*N3+xright*N2*N3+is*lsize);
                        x_r[rank] = xi[naxright];
                        for (bca=0; bca<N1/2; bca++) {							//broadcast x_r[rank] from #bca to all processors
                            MPI_Bcast(&x_r[bca], 1, MPI_DOUBLE, bca, MPI_COMM_WORLD);
                        }
                        MPI_Barrier(MPI_COMM_WORLD);
                        x_l[0] = x_r[(N1/2-1)];
                        
                        for (bca=0; bca<(N1/2-1); bca++) {
                            x_l[bca+1] = x_r[bca];
                        }
                        
                        for (bca=0; bca<N1/2; bca++) {
                            rho[is] += sqrt((x_r[bca]-x_l[bca])*(x_r[bca]-x_l[bca])/(dx*dx)+(y_r[bca]-y_l[bca])*(y_r[bca]-y_l[bca])/(dx*dx))/(burger*burger)/(N1*N2*N3);
                        }
                        
                    }
                    
                    else if (i==(lN1-1)) {
                        xleft = i-dx/2;
                        naxleft = 2*(k+j*N3+xleft*N2*N3+is*lsize);
                        x_l[rank] = xi[naxleft];
                        for (bca=0; bca<N1/2; bca++) {							//broadcast x_l[rank] from #bca to all processors
                            MPI_Bcast(&x_l[bca], 1, MPI_DOUBLE, bca, MPI_COMM_WORLD);
                        }
                        MPI_Barrier(MPI_COMM_WORLD);
                        x_r[(N1/2-1)] = x_l[0];
                        for (bca=0; bca<(N1/2-1); bca++) {
                            x_r[bca] = x_l[bca+1];
                        }
                        for (bca=0; bca<N1/2; bca++) {
                            rho[is] += sqrt((x_r[bca]-x_l[bca])*(x_r[bca]-x_l[bca])/(dx*dx)+(y_r[bca]-y_l[bca])*(y_r[bca]-y_l[bca])/(dx*dx))/(burger*burger)/(N1*N2*N3);
                        }
                    }
                    
                    else {
                        xright = i+dx/2;
                        xleft = i-dx/2;
                        naxright = 2*(k+j*N3+xright*N2*N3+is*lsize);
                        naxleft = 2*(k+j*N3+xleft*N2*N3+is*lsize);
                        x_r[rank] = xi[naxright];
                        for (bca=0; bca<N1/2; bca++) {							//broadcast x_r[rank] from #bca to all processors
                            MPI_Bcast(&x_r[bca], 1, MPI_DOUBLE, bca, MPI_COMM_WORLD);
                        }
                        MPI_Barrier(MPI_COMM_WORLD);
                        
                        x_l[rank] = xi[naxleft];
                        for (bca=0; bca<N1/2; bca++) {							//broadcast x_l[rank] from #bca to all processors
                            MPI_Bcast(&x_l[bca], 1, MPI_DOUBLE, bca, MPI_COMM_WORLD);
                        }
                        MPI_Barrier(MPI_COMM_WORLD);
                        
                        for (bca=0; bca<N1/2; bca++) {
                            rho[is] += sqrt((x_r[bca]-x_l[bca])*(x_r[bca]-x_l[bca])/(dx*dx)+(y_r[bca]-y_l[bca])*(y_r[bca]-y_l[bca])/(dx*dx))/(burger*burger)/(N1*N2*N3);
                        }
                    }
                }//end i j k
    }//end is
    
    rho_tot = rho[0] + rho[1] + rho[2];
    if (rank == 0) {
        printf("dislocation density %e %e %e %e\n", rho[0], rho[1], rho[2], rho_tot);
    }
    
    MPI_Barrier(MPI_COMM_WORLD);
    
    return;
}


/******************************************************************/
void avestrain(double **avepsd, double **avepst, double ***eps, double *xi, int nsize, int lsize, double **sigma, double S11, double S12, double S44, double mu, int lN1, double **ave_epsd, int rank, double **avepsts, int N1, int N2, int N3, int ND, int NS)
{
    
#define 	DELTA(i, j)   ((i==j)?1:0)
    
    int i, j, k, l, is, nb, ida, idb;
    double S[ND][ND][ND][ND];
    
    /*set S matrix*/
    for (i=0; i<ND; i++)
        for (j=0; j<ND; j++)
            for (k=0; k<ND; k++)
                for (l=0; l<ND; l++)
                {
                    S[i][j][k][l] = S44/4.0 * (DELTA(i,k)*DELTA(j,l)+DELTA(i,l)*DELTA(j,k))+S12*DELTA(i,j)*DELTA(k,l);
                }
    
    /*calculating average stress free strain, avepsd*/
    
    //for(is=0;is<NSV;is++)
    for(is=0;is<NS;is++)
    {
        for (ida=0; ida<ND; ida++)
            for (idb=0; idb<ND; idb++)
            {
                for(i=0;i<lN1;i++)
                    for(j=0;j<N2;j++)
                        for(k=0;k<N3;k++)
                        {
                            nb = 2*(k + j*N3 + i*N3*N2+ is*lsize);
                            avepsd[ida][idb]  += eps[is][ida][idb] * xi[nb];
                        }
                //printf("avepsd[ida][idb] %lf, rank %d, ida %d, idb %d\n", avepsd[ida][idb], rank, ida, idb);
                
                MPI_Allreduce(&avepsd[ida][idb], &ave_epsd[ida][idb], 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                
                //printf("ave_epsd[ida][idb] %lf, rank %d, ida %d, idb %d, ave_eps/nisze %lf\n", ave_epsd[ida][idb], rank, ida, idb, ave_epsd[ida][idb]/nsize);
                
                ave_epsd[ida][idb] = ave_epsd[ida][idb]/nsize;
            }
    }
    
    /*calculating microscopic strain, avepst*/
    
    for(i=0;i<ND;i++)
        for(j=0;j<ND;j++)
        {
            for(k=0;k<ND;k++)
                for(l=0;l<ND;l++)
                {
                    avepst[i][j] += S[i][j][k][l]*sigma[k][l]*mu;
                    avepsts[i][j] += S[i][j][k][l]*sigma[k][l]*mu;
                }
            avepst[i][j] += ave_epsd[i][j];
        }
    
    return;
}

/******************************************************************/

void strain(fftw_complex *data_strain, fftw_complex *data_disp, double *data_eps, fftw_complex *data_fftw, double *FF, double *UU, int tmstep, double **avepst, int lxs, int lN1, int nsize, int lsize, fftw_complex *work_strain, int rank)
{
    int i, j, k, l, is, na0, na1, nb, psys, ida, idb, index, index1, index2, nb1, nb2, nb3;
    int na11, na12, na13, na21, na22, na23, na31, na32, na33, ia, ib;
    double pi = 3.141592654;
    fftwnd_mpi_plan iiplan;
    FILE *of6;
    char outfile6[100], output6[100];
    
    
    iiplan = fftw3d_mpi_create_plan(MPI_COMM_WORLD, N1, N2, N3, FFTW_BACKWARD, FFTW_ESTIMATE);
    
    /*open output file for strain field*/
    
    for(is=0;is<NS;is++){
        strcpy(outfile6,scratch);
        ia = sprintf(output6, "/outstrain_it%08.0d_P%03.0d.dat", tmstep, rank);
        strcat(outfile6,output6);
        for(ib=0; ib<strlen(outfile6); ib++){
            if(outfile6[ib]==' ') outfile6[ib]='0';
        }
        
        of6 = fopen(outfile6,"w");
        
        if(rank==0){
            fprintf(of6,"zone  I = %d J = %d K = %d \n", N1, N2, N3);
        }
    }
    
    /*initialize*/
    
    for (i=0; i<lsize*ND*ND; i++){
        data_strain[i].re = 0.0;
        data_strain[i].im = 0.0;
    }
    for(i=0;i<lsize*ND; i++){
        data_disp[i].re = 0.0;
        data_disp[i].im = 0.0;
    }
    for (i=0; i<2*lsize*ND*ND; i++){
        data_eps[i] = 0.0;
    }
    
    /*calculate the total strain */
    
    for(is=0;is<NS;is++){
        for (ida=0; ida<ND; ida++)
            for (idb=0; idb<ND; idb++){
                for(i=0;i<lN1;i++)
                    for(j=0;j<N2;j++)
                        for(k=0;k<N3;k++){
                            index = k + j*N3 + i*N2*N3;
                            index1 = index + is*lsize;
                            index2 = index + ida*lsize + idb*lsize*ND;
                            nb = index1 + ida*lsize*NS + idb*lsize*NS*ND;
                            data_strain[index2].re += data_fftw[index1].re * FF[nb];
                            data_strain[index2].im += data_fftw[index1].im * FF[nb];
                        }
            }
    }
    
    for(is=0;is<NS;is++){
        for (ida=0; ida<ND; ida++){
            for(i=0;i<lN1;i++)
                for(j=0;j<N2;j++)
                    for(k=0;k<N3;k++){
                        index = k + j*N3 + i*N2*N3;
                        index1 = index + is*lsize;
                        index2 = index + ida*lsize;
                        nb = index1 + ida*lsize*NS;
                        data_disp[index2].re -= data_fftw[index1].im * UU[nb];	//disp=iUU/2pi
                        data_disp[index2].im += data_fftw[index1].re * UU[nb];
                    }
        }
    }
    
    
    for(i=0;i<ND;i++){
        psys = i*lsize;
        fftwnd_mpi(iiplan, 1, data_disp+psys, work_strain, FFTW_NORMAL_ORDER); // Inverse FFT (multiple)
        
    }
    for (i=0; i<lsize*ND; i++)
    {
        data_disp[i].re = data_disp[i].re/(nsize)/2/pi;
        data_disp[i].im = data_disp[i].im/(nsize)/2/pi;
    }
    
    for(i=0;i<ND;i++){
        for (j=0;j<ND;j++){
            psys = i*lsize + j*lsize*ND;
            fftwnd_mpi(iiplan, 1, data_strain+psys, work_strain, FFTW_NORMAL_ORDER); // Inverse FFT (multiple)
        }
    }
    
    //add in other two terms in strain
    
    for(ida=0;ida<ND;ida++)
        for (idb=0;idb<ND;idb++){
            for(i=0;i<lN1;i++)
                for(j=0;j<N2;j++)
                    for (k=0;k<N3;k++){
                        index = k + j*N3 + i*N2*N3 + ida*lsize + idb*lsize*ND;
                        if(load == 0) //stress control
                            data_strain[index].re = data_strain[index].re/nsize + avepst[ida][idb];
                        else if(load == 1) //strain control
                            data_strain[index].re = data_strain[index].re/nsize + epsilon[ida][idb];
                        
                        data_strain[index].im /= nsize;
                    }
        }
    
    for(ida=0;ida<ND;ida++)
        for (idb=0;idb<ND;idb++){
            for(i=0;i<lN1;i++)
                for(j=0;j<N2;j++)
                    for (k=0;k<N3;k++){
                        index = k + j*N3 + i*N2*N3;
                        index1 = index + ida*lsize + idb*lsize*ND;
                        index2 = index + idb*lsize + ida*lsize*ND;
                        na0 = 2*index1;
                        na1=na0+1;
                        data_eps[na0] = (data_strain[index1].re + data_strain[index2].re)/2.0;
                        data_eps[na1] = (data_strain[index1].im + data_strain[index2].im)/2.0;
                    }//end ijk
        }	 // end ida idb
    
    //if(it == NT-1 && (itp == (NSI -1))){
    for(i=0;i<lN1;i++)
        for(j=0;j<N2;j++)
            for(k=0;k<N3;k++){
                index = k + j*N3 + i*N2*N3;
                na0 = 2*index;  // index+i*N1*N2*N3+j*N1*N2*N3*ND);
                na11 = na0;
                na12 = na11 + 2*lsize;  //[1][2]
                na13 = na12 + 2*lsize;
                na21 = na13 + 2*lsize; //[2][1]
                na22 = na21 + 2*lsize;
                na23 = na22 + 2*lsize;
                na31 = na23 + 2*lsize;
                na32 = na31 + 2*lsize;
                na33 = na32 + 2*lsize;
                nb1 = index;
                nb2 = index + 1*lsize;
                nb3 = index + 2*lsize;
                
                fprintf(of6,"%d %d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %e %e %e  %e %e %e \n", lxs+i, j, k, data_eps[na11], data_eps[na12], data_eps[na13], data_eps[na21], data_eps[na22], data_eps[na23],  data_eps[na31], data_eps[na32], data_eps[na33], data_disp[nb1].re, data_disp[nb1].im, data_disp[nb2].re, data_disp[nb2].im, data_disp[nb3].re, data_disp[nb3].im);
                
            }
    
    fclose(of6);
    if(rank == 0){
        printf("leaving strain and displacement \n");
    }
    
    fftwnd_mpi_destroy_plan(iiplan);
    
    return;
}

/******************************************************************/

void stress(double *data_epsd, double *data_sigma, fftw_complex *data_strain, double *data_eps, double *xi, double ***eps, double C11, double C12, double C44, int tmstep, double **avesigma, double b, int lxs, int lN1, int nsize, int lsize, int rank, double **ave_sigma, double **sigma, double **avepsts, double *E_elas_test)
{
#define DELTA(i, j)   ((i==j)?1:0)
    
    int i, j, k, l, m, ida, idb, na, nb, na0, is, ia, ib, index;
    int na11, na12, na13, na21, na22, na23, na31, na32, na33;
    double C[ND][ND][ND][ND];
    double mu, xnu, young, ll;
    FILE *of7;
    char outfile7[100], output7[100];
    
    /* set Cijkl*/
    
    mu = C44-(2.0*C44+C12-C11)/5.0;
    ll = C12-(2.0*C44+C12-C11)/5.0;
    young = mu*(3*ll+2*mu)/(ll+mu);
    xnu = young/2.0/mu-1.0;
    
    /*open output file*/
    
    strcpy(outfile7,scratch);
    ia = sprintf(output7, "/outstress_it%08.0d_P%03.0d.dat", tmstep, rank);
    strcat(outfile7,output7);
    for(ib=0; ib<strlen(outfile7); ib++){
        if(outfile7[ib]==' ') outfile7[ib]='0';
    }
    
    of7 = fopen(outfile7,"w");
    
    if(rank==0){
        fprintf(of7,"zone  I = %d J = %d K = %d \n", N1, N2, N3);
    }
    
    for (i=0; i<ND; i++)
        for (j=0; j<ND; j++)
            for (k=0; k<ND; k++)
                for (m=0; m<ND; m++)
                {
                    C[i][j][k][m] = mu * (DELTA(i,k)*DELTA(j,m)+DELTA(i,m)*DELTA(j,k))+ll*DELTA(i,j)*DELTA(k,m);
                }
    
    for (i=0; i<2*lsize*ND*ND; i++)
    {
        data_sigma[i] =0;
        data_epsd[i]=0;
    }
    
    //for(is=0;is<NSV;is++)
    for(is=0;is<NS;is++)
    {
        for (ida=0; ida<ND; ida++)
            for (idb=0; idb<ND; idb++)
            {
                for(i=0;i<lN1;i++)
                    for(j=0;j<N2;j++)
                        for(k=0;k<N3;k++)
                        {
                            na = 2*(i*N2*N3 + j*N3 + k + ida*lsize + idb*lsize*ND);
                            nb = 2*(k + j*N3 + i*N3*N2 + is*lsize);
                            data_epsd[na] += eps[is][ida][idb] * xi[nb];
                        }
            }
    }
    
    *E_elas_test = 0.0;
    for (ida=0;ida<ND;ida++)
        for (idb=0;idb<ND;idb++)
        {
            for(i=0;i<lN1;i++)
                for(j=0;j<N2;j++)
                    for(k=0;k<N3;k++)
                    {
                        index = k + j*N3 + i*N2*N3;
                        na = 2*(index + ida*lsize + idb*lsize*ND);
                        for (m=0;m<ND;m++)
                            for (l=0;l<ND;l++)
                            {
                                nb = 2*(index + m*lsize + l*lsize*ND);
                                
                                if(load == 0)
                                    data_sigma[na] += C[ida][idb][m][l]*(data_eps[nb] - data_epsd[nb] - avepsts[m][l]);
                                else {
                                    data_sigma[na] += C[ida][idb][m][l]*(data_eps[nb] - data_epsd[nb]);
                                }
                                *E_elas_test += C[ida][idb][m][l]*((data_eps[na] - data_epsd[na])*(data_eps[nb] - data_epsd[nb]) + (data_eps[na+1] - data_epsd[na+1])*(data_eps[nb+1] - data_epsd[nb+1]));
                                
                                data_sigma[na+1] = 0.0;
                            }
                        if(load == 0)
                            data_sigma[na]+=sigma[ida][idb]*mu;
                        avesigma[ida][idb] += data_sigma[na];
                    }
            
            MPI_Reduce(&avesigma[ida][idb], &ave_sigma[ida][idb], 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
            
            if(rank == 0){
                ave_sigma[ida][idb] = ave_sigma[ida][idb]/nsize;
            }
        } // end ida idb
    MPI_Allreduce(E_elas_test, E_elas_test, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    *E_elas_test *= 0.5*b*b/N1;
    
    //if(it == NT-1 && (itp == (NSI -1))){
    for(i=0;i<lN1;i++)
        for(j=0;j<N2;j++)
            for(k=0;k<N3;k++){
                na0 = 2*(k + j*N3 + i*N2*N3);  // +i*N1*N2*N3+j*N1*N2*N3*ND);
                na11 = na0;
                na12 = na11 + 2*lsize;  //[1][2]
                na13 = na12 + 2*lsize;
                na21 = na13 + 2*lsize; //[2][1]
                na22 = na21 + 2*lsize;
                na23 = na22 + 2*lsize;
                na31 = na23 + 2*lsize;
                na32 = na31 + 2*lsize;
                na33 = na32 + 2*lsize;
                fprintf(of7,"%d %d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf\n", lxs+i, j, k, data_sigma[na11], data_sigma[na12], data_sigma[na13], data_sigma[na21], data_sigma[na22], data_sigma[na23], data_sigma[na31], data_sigma[na32], data_sigma[na33]);
            }
    
    if(rank == 0){
        printf("E_elas_test %e leaving stress\n", *E_elas_test);
    }
    fclose(of7);
    return;
}
/******************************************************************/

void core_ener(int cflag, int it, int rank, int lxs, double An, int lN1, double c0, double c1, double c2, double c3, double c4, double a1, double a3, double **xn, double **xb, fftw_complex *data_fftw, fftw_complex *data_core, double dslip, double b, double *fcore, double *df1core, double *df2core, double *df3core, double *dE_core, double E_core, int itp, double pi, int N1, int N2, int N3, int ND, int NS, int NP, int NT, int NSI, double Bn, double mu, int nsize, int lsize, double size, double isf, char *scratch,double usf)
{
    /*cflag == 1 -> assume perfect dislocations only, none extended
     cflag == 2 -> models all dislocations as extended, incorporates partials
     cflag == 3 -> extended dislocations sine approximation using intrinsic and unstable stacking fault energies as only parameters -- 1D version
     cflag == 4 -> extended dislocations approximated with only isf and usf--3D version*/
    
    int i, j, k, isa, index, index1, index2, index3, indexdx, plane, num, ia, ib, ic, id, ie, ig, ih, ij;
    int counter, marker, indexmin, tag;
    double *delta, *ddelta, dx, mpidel, *f_core, p;
    FILE *of4, *of5 /**of6, *of7*/;
    char outfile4[100], output4[100], outfile5[100], output5[100]/*outfile6[100], output6[100], outfile7[100], output7[100]*/;
    MPI_Status status;
    
    delta = (double*) malloc((NP)*lsize*sizeof(double));
    //ddelta = (double*) malloc((NP)*lsize*sizeof(double));
    ddelta = (double*) malloc((NP)*(lN1)*sizeof(double));
    f_core = (double*) malloc((NP)*lsize*sizeof(double));
    
    dx = size/N1;
    tag = 1;
    c0 = c0/(mu);
    c1 = c1/(mu);
    c2 = c2/(mu);
    c3 = c3/(mu);
    c4 = c4/(mu);
    a1 = a1/(mu);
    a3 = a3/(mu);
    isf = isf/(mu*dslip*b);
    usf = usf/(mu*dslip*b);
    An = An/(mu*dslip*b);
    p = pi/(sqrt(3.0)*(b/1.0E-10));
    
    //mark: define a different intrinsic stacking fault and unstable stacking fault for cflag == 3
    double isf_local, usf_local;
    isf_local = 84.718E-3/(mu*dslip*b);
    usf_local = 211.688E-3/(mu*dslip*b);
    
    /*if(rank == 0){
     strcpy(outfile7,scratch);
     ih = sprintf(output7, "/outcore_it%05.0d.dat", it);
     strcat(outfile7,output7);
     
     for(ij=0; ij<strlen(outfile7); ij++){
     if(outfile7[ij]==' ') outfile7[ij]='0';
     }
     
     of7 = fopen(outfile7,"w");
     }*/
    
    //if(fmod((double)(it),1000.0) == 0.0 || it == NT-1){
    /*if((fmod((double)(itp),5.0) == 0.0 && it == NT-1) || it == NT-1 && (itp == (NSI -1))){
     strcpy(outfile4,scratch);
     ia = sprintf(output4, "/outdelta_it%05.0d_P%03.0d.dat",it, rank);
     strcat(outfile4,output4);
     
     for(ib=0; ib<strlen(outfile4); ib++){
     if(outfile4[ib]==' ') outfile4[ib]='0';
     }
     
     of4 = fopen(outfile4,"w");
     
     strcpy(outfile5,scratch);
     ic = sprintf(output5, "/outdeltadx_it%05.0d_P%03.0d.dat", it, rank);
     strcat(outfile5,output5);
     
     for(id=0; id<strlen(outfile5); id++){
     if(outfile5[id]==' ') outfile5[id]='0';
     }
     
     of5 = fopen(outfile5,"w");
     
     if(rank == 0){
     fprintf(of4,"zone   I = %d\n", N1);
     fprintf(of5,"zone   I = %d\n", N1);
     }
     }/*it*/
    
    if(cflag == 1){
        for(isa=0;isa<NS;isa++)
        {
            for(i=0;i<lN1;i++)
                for(j=0;j<N2;j++)
                    for(k=0;k<N3;k++)
                    {
                        index = i*N2*N3 + j*N3 + k + isa*lsize;
                        E_core += An*(sin(pi*data_fftw[index].re)*sin(pi*data_fftw[index].re))/N1;
                        dE_core[index] = An*pi*sin(2.0*pi*data_fftw[index].re);
                        //if(rank == 0){
                        //printf("%lf, %lf, %lf\n", E_core, data_fftw[index].re, pi);
                        //}
                        //printf("%d %d %d %d  %lf %lf\n",i, j, k, ,index, dE_core[index], E_core);
                    }/*ijk*/
            //if(rank == 0){
            //printf("%d %d %d %d %lf %lf\n",i, j, k, index, dE_core[index], E_core);
            //}
        }/*isa*/
        if(rank == 0 && it == NT-1){
            printf("Core Energy for this time step %lf\n", E_core);
        }
        
        /*if(it == NT-1 && (itp == (NSI -1))){
         for(i=0;i<lsize*NP;i++){
         delta[i] = 0.0;
         }
         for(plane=0;plane<NP;plane++)
         {
         for(i=0;i<lN1;i++)
         for(j=0;j<N2;j++)
         for(k=0;k<N3;k++)
         {
         index = i*N2*N3 + j*N3 + k + plane*lsize;
         //indexdx = ((lxs+i)+1)*N2*N3 + j*N3 + k + plane*lsize;
         index1 = i*N2*N3 + j*N3 + k + 0*lsize + plane*lsize*3;
         index2 = i*N2*N3 + j*N3 + k + 1*lsize + plane*lsize*3;
         index3 = i*N2*N3 + j*N3 + k + 2*lsize + plane*lsize*3;
         
         //delta[index] = (data_fftw[index1].re*xb[0][0] + data_fftw[index2].re*xb[1][0] + data_fftw[index3].re*xb[2][0])*xb[1][0] + (data_fftw[index1].re*xb[0][1] + data_fftw[index2].re*xb[1][1] + data_fftw[index3].re*xb[2][1])*xb[1][1] + (data_fftw[index1].re*xb[0][2] + data_fftw[index2].re*xb[1][2] + data_fftw[index3].re*xb[2][2])*xb[1][2];
         
         //if(lxs+i == N1/2 && k == N3/2){
         if(j == N2/2 && k == N3/2){
         //if(rank == 0){
         //fprintf(of4,"zone   J = %d\n", N2);
         //}
         fprintf(of4, "%d %lf %lf\n", lxs+i, (dx)*(lxs+i), delta[index]);
         }
         }//ijk
         }//plane
         
         for(plane=0;plane<NP;plane++)
         {
         for(i=0;i<lN1;i++)
         //for(j=0;j<N2;j++)
         //for(k=0;k<N3;k++)
         {
         //if(j == N2/2 && k == N3/2){
         j = N2/2;
         k = N3/2;
         index = i*N2*N3 + j*N3 + k + plane*lsize;
         //indexdx = i*N2*N3 + (j+1)*N3 + k + plane*lsize;
         indexdx = (i+1)*N2*N3 + j*N3 + k + plane*lsize;
         if((i+1) == lN1){
         indexmin = 0*N2*N3 + j*N3 + k + plane*lsize;
         if(rank != 0){
         MPI_Send(&delta[indexmin], 1, MPI_DOUBLE, rank-1, tag, MPI_COMM_WORLD);
         }
         if((lxs+(i+1)) == N1){
         mpidel = delta[index];
         ddelta[i] = (mpidel - delta[index])/dx;
         goto skip_point;
         }
         MPI_Recv(&mpidel, 1, MPI_DOUBLE, rank+1, tag, MPI_COMM_WORLD, &status);
         ddelta[i] = (mpidel - delta[index])/dx;
         }
         else{
         ddelta[i] = (delta[indexdx] - delta[index])/dx;
         }
         
         //ddelta[index] = (delta[indexdx] - delta[index])/dx;
         //if(lxs+i == N1/2 && k == N3/2){
         //if(rank == 0){
         //fprintf(of5,"zone   J = %d", N2);
         //}
         //fprintf(of5, "%d %lf %lf %d %d %lf %lf %lf\n", j, (dslip/N1)*j, ddelta[index], index, indexdx, delta[index], delta[indexdx], dx);
         
         skip_point:
         continue;
         //fprintf(of5, "%d %lf %lf\n", lxs+i, (dx)*(lxs+i), ddelta[i]);
         }
         }
         }//it*/
        
    }//cflag
    
    else if(cflag == 2){
        
        counter = 0;
        marker = 0;
        
        /*as of now the planes need to be in the correct order initially for this to work
         also make sure to call this subroutine before the FFTW functions are called so that data_fftw = xi[na0] or double check that this is true.
         subroutine calculates core energy for each time step, need to called every time step and initialize every time step.*/
        
        for(isa=0;isa<NS;isa++)
        {
            for(i=0;i<lN1;i++)
                for(j=0;j<N2;j++)
                    for(k=0;k<N3;k++)
                    {
                        index = i*N2*N3 + j*N3 + k + isa*lsize;
                        data_core[index].re = data_fftw[index].re;
                        data_core[index].im = data_fftw[index].im;
                        //if((lxs+i) >= N1/2 && k == N3/2){
                        //printf("%d %d %d %d core %lf fftw %lf index %d\n", i, j, k, isa, data_core[index].re, data_fftw[index].re, index);
                        //}
                    }
        }
        
        for(plane=0;plane<NP;plane++)
        {
            for(i=0;i<lN1;i++)
                for(j=0;j<N2;j++)
                    for(k=0;k<N3;k++)
                    {
                        index = i*N2*N3 + j*N3 + k + plane*lsize;
                        index1 = i*N2*N3 + j*N3 + k + 0*lsize + plane*lsize*3;
                        index2 = i*N2*N3 + j*N3 + k + 1*lsize + plane*lsize*3;
                        index3 = i*N2*N3 + j*N3 + k + 2*lsize + plane*lsize*3;
                        
                        /*no derivatives taken use to calculate E_core*/
                        
                        fcore[index] = (c0 + c1*(cos(2.0*pi*(data_core[index1].re-data_core[index2].re)) + cos(2.0*pi*(data_core[index2].re-data_core[index3].re)) + cos(2.0*pi*(data_core[index3].re-data_core[index1].re))) + c2*(cos(2.0*pi*(2.0*data_core[index1].re-data_core[index2].re-data_core[index3].re)) + cos(2.0*pi*(2.0*data_core[index2].re-data_core[index3].re-data_core[index1].re)) + cos(2.0*pi*(2.0*data_core[index3].re-data_core[index1].re-data_core[index2].re))) + c3*(cos(4.0*pi*(data_core[index1].re-data_core[index2].re)) + cos(4.0*pi*(data_core[index2].re-data_core[index3].re)) + cos(4.0*pi*(data_core[index3].re-data_core[index1].re))) + c4*(cos(2.0*pi*(3.0*data_core[index1].re-data_core[index2].re-2.0*data_core[index3].re)) + cos(2.0*pi*(3.0*data_core[index1].re-2.0*data_core[index2].re-data_core[index3].re)) + cos(2*pi*(3.0*data_core[index2].re-data_core[index3].re-2.0*data_core[index1].re)) + cos(2.0*pi*(3.0*data_core[index2].re-2.0*data_core[index3].re-data_core[index1].re)) + cos(2.0*pi*(3.0*data_core[index3].re-data_core[index1].re-2.0*data_core[index2].re)) + cos(2.0*pi*(3.0*data_core[index3].re-2.0*data_core[index1].re-data_core[index2].re))) + a1*(sin(2.0*pi*(data_core[index1].re-data_core[index2].re)) + sin(2.0*pi*(data_core[index2].re-data_core[index3].re)) + sin(2.0*pi*(data_core[index3].re-data_core[index1].re))) + a3*(sin(4.0*pi*(data_core[index1].re-data_core[index2].re)) + sin(4.0*pi*(data_core[index2].re-data_core[index3].re)) + sin(4.0*pi*(data_core[index3].re-data_core[index1].re))))/(dslip*b);
                        
                        //MPI_Reduce(&fcore[index], &f_core[index], 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
                        
                        E_core += fcore[index];
                        //if(rank == 0){
                        //  E_core += Bn*(f_core[index]/nsize);
                        //}
                        
                        /*partial derivative wrt phase field 1*/
                        
                        
                        df1core[index] = ((2.0*pi)/(dslip*b))*(c1*(sin(2.0*pi*(data_core[index2].re-data_core[index1].re)) + sin(2.0*pi*(data_core[index3].re-data_core[index1].re))) + c2*(2.0*sin(2.0*pi*(data_core[index2].re+data_core[index3].re-2.0*data_core[index1].re)) + sin(2.0*pi*(2.0*data_core[index2].re-data_core[index3].re-data_core[index1].re)) + sin(2.0*pi*(2.0*data_core[index3].re-data_core[index1].re-data_core[index2].re))) + (2.0*c3)*(sin(4.0*pi*(data_core[index2].re-data_core[index1].re)) + sin(4.0*pi*(data_core[index3].re-data_core[index1].re))) + c4*(3.0*sin(2.0*pi*(data_core[index2].re+2.0*data_core[index3].re-3.0*data_core[index1].re)) + 3.0*sin(2.0*pi*(2.0*data_core[index2].re+data_core[index3].re-3.0*data_core[index1].re)) + 2.0*sin(2.0*pi*(3.0*data_core[index2].re-data_core[index3].re-2.0*data_core[index1].re)) + sin(2.0*pi*(3.0*data_core[index2].re-2.0*data_core[index3].re-data_core[index1].re)) + sin(2.0*pi*(3.0*data_core[index3].re-data_core[index1].re-2.0*data_core[index2].re)) + 2.0*sin(2.0*pi*(3.0*data_core[index3].re-2.0*data_core[index1].re-data_core[index2].re))) + a1*(cos(2.0*pi*(data_core[index1].re-data_core[index2].re)) - cos(2.0*pi*(data_core[index3].re-data_core[index1].re))) + (2.0*a3)*(cos(4.0*pi*(data_core[index1].re-data_core[index2].re)) - cos(4.0*pi*(data_core[index3].re-data_core[index1].re))));
                        
                        /*partial derivative wrt phase field 2*/
                        
                        
                        df2core[index] = ((2.0*pi)/(dslip*b))*(c1*(sin(2.0*pi*(data_core[index1].re-data_core[index2].re)) + sin(2.0*pi*(data_core[index3].re-data_core[index2].re))) + c2*(sin(2.0*pi*(2.0*data_core[index1].re-data_core[index2].re-data_core[index3].re)) + 2.0*sin(2.0*pi*(data_core[index3].re+data_core[index1].re-2.0*data_core[index2].re)) + sin(2.0*pi*(2.0*data_core[index3].re-data_core[index1].re-data_core[index2].re))) + (2.0*c3)*(sin(4.0*pi*(data_core[index1].re-data_core[index2].re)) + sin(4.0*pi*(data_core[index3].re-data_core[index2].re))) + c4*(sin(2.0*pi*(3.0*data_core[index1].re-data_core[index2].re-2.0*data_core[index3].re)) + 2.0*sin(2.0*pi*(3.0*data_core[index1].re-2.0*data_core[index2].re-data_core[index3].re)) + 3.0*sin(2.0*pi*(data_core[index3].re+2.0*data_core[index1].re-3.0*data_core[index2].re)) + 3.0*sin(2.0*pi*(2.0*data_core[index3].re+data_core[index1].re-3.0*data_core[index2].re)) + 2.0*sin(2.0*pi*(3.0*data_core[index3].re-data_core[index1].re-2.0*data_core[index2].re)) + sin(2.0*pi*(3.0*data_core[index3].re-2.0*data_core[index1].re-data_core[index2].re))) + a1*(cos(2.0*pi*(data_core[index2].re-data_core[index3].re)) - cos(2.0*pi*(data_core[index1].re-data_core[index2].re))) + (2.0*a3)*(cos(4.0*pi*(data_core[index2].re-data_core[index3].re)) - cos(4.0*pi*(data_core[index1].re-data_core[index2].re))));
                        
                        /*partial derivative wrt phase field 3*/
                        
                        
                        df3core[index] = ((2.0*pi)/(dslip*b))*(c1*(sin(2.0*pi*(data_core[index2].re-data_core[index3].re)) + sin(2.0*pi*(data_core[index1].re-data_core[index3].re))) + c2*(sin(2.0*pi*(2.0*data_core[index1].re-data_core[index2].re-data_core[index3].re)) + sin(2.0*pi*(2.0*data_core[index2].re-data_core[index3].re-data_core[index1].re)) + 2.0*sin(2.0*pi*(data_core[index1].re+data_core[index2].re-2.0*data_core[index3].re))) + (2.0*c3)*(sin(4.0*pi*(data_core[index2].re-data_core[index3].re)) + sin(4.0*pi*(data_core[index1].re-data_core[index3].re))) + c4*(2.0*sin(2.0*pi*(3.0*data_core[index1].re-data_core[index2].re-2.0*data_core[index3].re)) + sin(2.0*pi*(3.0*data_core[index1].re-2.0*data_core[index2].re-data_core[index3].re)) + sin(2.0*pi*(3.0*data_core[index2].re-data_core[index3].re-2.0*data_core[index1].re)) + 2.0*sin(2.0*pi*(3.0*data_core[index2].re-2.0*data_core[index3].re-data_core[index1].re)) + 3.0*sin(2.0*pi*(data_core[index1].re+2.0*data_core[index2].re-3.0*data_core[index3].re)) + 3.0*sin(2.0*pi*(2.0*data_core[index1].re+data_core[index2].re-3.0*data_core[index3].re))) + a1*(cos(2.0*pi*(data_core[index3].re-data_core[index1].re)) - cos(2.0*pi*(data_core[index2].re-data_core[index3].re))) + (2.0*a3)*(cos(4.0*pi*(data_core[index3].re-data_core[index1].re)) - cos(4.0*pi*(data_core[index2].re-data_core[index3].re))));
                        
                        dE_core[index1] = Bn*df1core[index];
                        dE_core[index2] = Bn*df2core[index];
                        dE_core[index3] = Bn*df3core[index];
                        
                    }//ijk
        }//plane
        
        MPI_Reduce(&E_core, &E_core, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if(rank == 0){
            E_core *= Bn;
        }
        //if(rank == 0){
        //printf("fcore %lf df1core %lf df2core %lf df3core %lf index %d\n", fcore[index], df1core[index], df2core[index], df3core[index], index);
        //}
        /*for(plane=0;plane<NP;plane++)
         {
         for(i=0;i<lN1;i++)
         for(j=0;j<N2;j++)
         for(k=0;k<N3;k++)
         {
         index = i*N2*N3 + j*N3 + k + plane*lsize;
         index1 = i*N2*N3 + j*N3 + k + 0*lsize + plane*lsize*3;
         index2 = i*N2*N3 + j*N3 + k + 1*lsize + plane*lsize*3;
         index3 = i*N2*N3 + j*N3 + k + 2*lsize + plane*lsize*3;
         
         dE_core[index1] = Bn*df1core[index];
         dE_core[index2] = Bn*df2core[index];
         dE_core[index3] = Bn*df3core[index];
         }
         }*/
        
        if(rank == 0 && it == NT-1){
            printf("Core Energy for this time step %lf\n", E_core);
            //fprintf(of7, "%d %lf\n", it, E_core);
        }
        /*if(fmod((double)(it),1000.0) == 0.0 || it == NT-1){
         for(isa=0;isa<NS;isa++)
         {
         strcpy(outfile6,scratch);
         ie = sprintf(output6, "/outcorevar_it%05.0d_NS%02.0d_P%03.0d.dat", it, isa, rank);
         strcat(outfile6,output6);
         
         for(ig=0; ig<strlen(outfile6); ig++){
         if(outfile6[ig]==' ') outfile6[ig]='0';
         }
         
         of6 = fopen(outfile6,"w");
         
         if(rank == 0){
         fprintf(of6,"zone   I = %d\n", N1);
         }
         
         for(i=0;i<lN1;i++)
         for(j=0;j<N2;j++)
         for(k=0;k<N3;k++)
         {
         index = i*N2*N3 + j*N3 + k + isa*lsize;
         if(j == N2/2 && k == N3/2){
         fprintf(of6, "%d %lf %lf\n", lxs+i, (dx)*(lxs+i), dE_core[index]);
         }
         } //ijk
         }
         }*/
        
        //MPI_Barrier(MPI_COMM_WORLD);
        //if(it == NT-1){
        //if(fmod((double)(it),1000.0) == 0.0 || it == NT-1){
        /*if((fmod((double)(itp),5.0) == 0.0 && it == NT-1) || it == NT-1 && (itp == (NSI -1))){
         for(i=0;i<lsize*NP;i++){
         delta[i] = 0.0;
         }
         for(plane=0;plane<NP;plane++)
         {
         for(i=0;i<lN1;i++)
         for(j=0;j<N2;j++)
         for(k=0;k<N3;k++)
         {
         index = i*N2*N3 + j*N3 + k + plane*lsize;
         //indexdx = ((lxs+i)+1)*N2*N3 + j*N3 + k + plane*lsize;
         index1 = i*N2*N3 + j*N3 + k + 0*lsize + plane*lsize*3;
         index2 = i*N2*N3 + j*N3 + k + 1*lsize + plane*lsize*3;
         index3 = i*N2*N3 + j*N3 + k + 2*lsize + plane*lsize*3;
         
         delta[index] = (data_core[index1].re*xb[0][0] + data_core[index2].re*xb[1][0] + data_core[index3].re*xb[2][0])*xb[1][0] + (data_core[index1].re*xb[0][1] + data_core[index2].re*xb[1][1] + data_core[index3].re*xb[2][1])*xb[1][1] + (data_core[index1].re*xb[0][2] + data_core[index2].re*xb[1][2] + data_core[index3].re*xb[2][2])*xb[1][2];
         
         //if(lxs+i == N1/2 && k == N3/2){
         if(j == N2/2 && k == N3/2){
         //if(rank == 0){
         //fprintf(of4,"zone   J = %d\n", N2);
         //}
         fprintf(of4, "%d %lf %lf\n", lxs+i, (dx)*(lxs+i), delta[index]);
         }
         } //ijk
         }//plane
         
         for(plane=0;plane<NP;plane++)
         {
         for(i=0;i<lN1;i++)
         //for(j=0;j<N2;j++)
         //for(k=0;k<N3;k++)
         {
         //if(j == N2/2 && k == N3/2){
         j = N2/2;
         k = N3/2;
         index = i*N2*N3 + j*N3 + k + plane*lsize;
         //indexdx = i*N2*N3 + (j+1)*N3 + k + plane*lsize;
         indexdx = (i+1)*N2*N3 + j*N3 + k + plane*lsize;
         if((i+1) == lN1){
         indexmin = 0*N2*N3 + j*N3 + k + plane*lsize;
         if(rank != 0){
         MPI_Send(&delta[indexmin], 1, MPI_DOUBLE, rank-1, tag, MPI_COMM_WORLD);
         }
         if((lxs+(i+1)) == N1){
         mpidel = delta[index];
         ddelta[i] = (mpidel - delta[index])/dx;
         goto skip_point1;
         }
         MPI_Recv(&mpidel, 1, MPI_DOUBLE, rank+1, tag, MPI_COMM_WORLD, &status);
         ddelta[i] = (mpidel - delta[index])/dx;
         }
         else{
         ddelta[i] = (delta[indexdx] - delta[index])/dx;
         }
         
         //ddelta[index] = (delta[indexdx] - delta[index])/dx;
         //if(lxs+i == N1/2 && k == N3/2){
         //if(rank == 0){
         //fprintf(of5,"zone   J = %d", N2);
         //}
         //fprintf(of5, "%d %lf %lf %d %d %lf %lf %lf\n", j, (dslip/N1)*j, ddelta[index], index, indexdx, delta[index], delta[indexdx], dx);
         
         skip_point1:
         fprintf(of5, "%d %lf %lf\n", lxs+i, (dx)*(lxs+i), ddelta[i]);
         }//ijk
         }//plane
         } //it    */
    }/*cflag*/
    
    else if(cflag == 3){
        for(isa=0;isa<NS;isa++)
        {
            for(i=0;i<lN1;i++)
                for(j=0;j<N2;j++)
                    for(k=0;k<N3;k++)
                    {
                        index = i*N2*N3 + j*N3 + k + isa*lsize;
                       // if (lxs+i<=N1*2/3.&&lxs+i>=N1*1/3.&&(j>=N2*3/4.||j<=N2*1/4.)) {
                         //if (((lxs+i)-N1*0.45)*((lxs+i)-N1*0.45)+(j-N2*0.)*(j-N2*0.) <= (N2*0.275)*(N2*0.275) || ((lxs+i)-N1*0.45)*((lxs+i)-N1*0.45)+(j-N2*1)*(j-N2*1) <= (N2*0.275)*(N2*0.275)) { //d=0.45N2
                        // if (((lxs+i)-N1*0.45)*((lxs+i)-N1*0.45)+(j-N2*0.)*(j-N2*0.) <= (N2*0.3)*(N2*0.3) || ((lxs+i)-N1*0.45)*((lxs+i)-N1*0.45)+(j-N2*1)*(j-N2*1) <= (N2*0.3)*(N2*0.3)) { //d=0.4N2
                        // if (((lxs+i)-N1*0.45)*((lxs+i)-N1*0.45)+(j-N2*0.)*(j-N2*0.) <= (N2*0.25)*(N2*0.25) || ((lxs+i)-N1*0.45)*((lxs+i)-N1*0.45)+(j-N2*1.)*(j-N2*1.) <= (N2*0.25)*(N2*0.25)) { //d=0.5N2
                        // if (((lxs+i)-N1*0.45)*((lxs+i)-N1*0.45)+(j-N2*0.)*(j-N2*0.) <= (N2*0.325)*(N2*0.325) || ((lxs+i)-N1*0.45)*((lxs+i)-N1*0.45)+(j-N2*1)*(j-N2*1) <= (N2*0.325)*(N2*0.325)) { //d=0.35N2
                        // if (((lxs+i)-N1*0.45)*((lxs+i)-N1*0.45)+(j-N2*0.)*(j-N2*0.) <= (N2*0.35)*(N2*0.35) || ((lxs+i)-N1*0.45)*((lxs+i)-N1*0.45)+(j-N2*1)*(j-N2*1) <= (N2*0.35)*(N2*0.35)) { //d=0.3N2
                          if (((lxs+i)-N1*0.5375)*((lxs+i)-N1*0.5375)+(j-N2*0.)*(j-N2*0.) <= (N2*0.425)*(N2*0.425) || ((lxs+i)-N1*0.5375)*((lxs+i)-N1*0.5375)+(j-N2*1)*(j-N2*1) <= (N2*0.425)*(N2*0.425)) { //d=0.3N2
 
                            isf_local = 0.5*isf;
                            usf_local = usf;

                            E_core+=(isf_local*(sin(pi*data_fftw[index].re)*sin(pi*data_fftw[index].re))+(usf_local-isf_local/2.0)*(sin(2*pi*data_fftw[index].re)*sin(2*pi*data_fftw[index].re)))/N1;
                            dE_core[index] = Bn*(isf_local*pi*sin(2*pi*data_fftw[index].re) + (usf_local-isf_local/2.0)*2*pi*sin(4*pi*data_fftw[index].re));
                        }else{
                            E_core += (isf*(sin(pi*data_fftw[index].re)*sin(pi*data_fftw[index].re)) + (usf-isf/2.0)*(sin(2*pi*data_fftw[index].re)*sin(2*pi*data_fftw[index].re)))/N1;
                            dE_core[index] = Bn*(isf*pi*sin(2*pi*data_fftw[index].re) + (usf-isf/2.0)*2*pi*sin(4*pi*data_fftw[index].re));
                        }
                    }/*ijk*/
        }/*isa*/
        if(rank == 0 && /*fmod((double)(it),1000.0) == 0.0*/ it == NT-1){
            printf("Core Energy for this time step %lf\n", E_core);
        }
        
        //if(fmod((double)(it),1000.0) == 0.0 || it == NT-1){
//        if(/*fmod((double)(itp),5.0) == 0.0 ||*/ it == NT-1 && (itp == (NSI -1))){
//            for(i=0;i<lsize*NP;i++){
//                delta[i] = 0.0;
//            }
//            for(plane=0;plane<NP;plane++)
 //           {
//                for(i=0;i<lN1;i++)
//                    for(j=0;j<N2;j++)
//                        for(k=0;k<N3;k++)
//                        {
//                            index = i*N2*N3 + j*N3 + k + plane*lsize;
//                            //indexdx = ((lxs+i)+1)*N2*N3 + j*N3 + k + plane*lsize;
//                            index1 = i*N2*N3 + j*N3 + k + 0*lsize + plane*lsize*3;
//                            index2 = i*N2*N3 + j*N3 + k + 1*lsize + plane*lsize*3;
//                            index3 = i*N2*N3 + j*N3 + k + 2*lsize + plane*lsize*3;
//
//                            delta[index] = /*(1.0/2.0)**/(data_fftw[index1].re*xb[0][0] + data_fftw[index2].re*xb[1][0] + data_fftw[index3].re*xb[2][0])*xb[1][0] + (data_fftw[index1].re*xb[0][1] + data_fftw[index2].re*xb[1][1] + data_fftw[index3].re*xb[2][1])*xb[1][1] + (data_fftw[index1].re*xb[0][2] + data_fftw[index2].re*xb[1][2] + data_fftw[index3].re*xb[2][2])*xb[1][2];
                            
                            //if(lxs+i == N1/2 && k == N3/2){
//                            if(j == N2/2 && k == N3/2){
                                //if(rank == 0){
                                //fprintf(of4,"zone   J = %d\n", N2);
                                //}
//                                fprintf(of4, "%d %lf %lf\n", lxs+i, (dx)*(lxs+i), delta[index]);
//                            }
//                        }/*ijk*/
//            }/*plane*/
            
//            for(plane=0;plane<NP;plane++)
//            {
//                for(i=0;i<lN1;i++)
                    //for(j=0;j<N2;j++)
                    //for(k=0;k<N3;k++)
//                {
                    //if(j == N2/2 && k == N3/2){
//                    j = N2/2;
//                    k = N3/2;
//                    index = i*N2*N3 + j*N3 + k + plane*lsize;
                    //indexdx = i*N2*N3 + (j+1)*N3 + k + plane*lsize;
//                    indexdx = (i+1)*N2*N3 + j*N3 + k + plane*lsize;
//                    if((i+1) == lN1){
//                        indexmin = 0*N2*N3 + j*N3 + k + plane*lsize;
//                        if(rank != 0){
//                            MPI_Send(&delta[indexmin], 1, MPI_DOUBLE, rank-1, tag, MPI_COMM_WORLD);
//                        }
//                        if((lxs+(i+1)) == N1){
//                            mpidel = delta[index];
//                            ddelta[i] = (mpidel - delta[index])/dx;
//                            goto skip_point2;
//                        }
//                        MPI_Recv(&mpidel, 1, MPI_DOUBLE, rank+1, tag, MPI_COMM_WORLD, &status);
//                        ddelta[i] = (mpidel - delta[index])/dx;
//                    }
//                    else{
//                        ddelta[i] = (delta[indexdx] - delta[index])/dx;
//                    }
                    
                    //ddelta[index] = (delta[indexdx] - delta[index])/dx;
                    //if(lxs+i == N1/2 && k == N3/2){
                    //if(rank == 0){
                    //fprintf(of5,"zone   J = %d", N2);
                    //}
                    //fprintf(of5, "%d %lf %lf %d %d %lf %lf %lf\n", j, (dslip/N1)*j, ddelta[index], index, indexdx, delta[index], delta[indexdx], dx);
                    
//                skip_point2:
//                    fprintf(of5, "%d %lf %lf\n", lxs+i, (dx)*(lxs+i), ddelta[i]);
//                }
//            }
//        }/*it*/
    }/*cflag*/
    
    else if(cflag == 4){ //not working 11/08/10
        for(isa=0;isa<NS;isa++)
        {
            for(i=0;i<lN1;i++)
                for(j=0;j<N2;j++)
                    for(k=0;k<N3;k++)
                    {
                        index = i*N2*N3 + j*N3 + k + isa*lsize;
                        data_core[index].re = data_fftw[index].re;
                        data_core[index].im = data_fftw[index].im;
                        //if((lxs+i) >= N1/2 && k == N3/2){
                        //printf("%d %d %d %d core %lf fftw %lf index %d\n", i, j, k, isa, data_core[index].re, data_fftw[index].re, index);
                        //}
                    }
        }
        
        for(plane=0;plane<NP;plane++)
        {
            for(i=0;i<lN1;i++)
                for(j=0;j<N2;j++)
                    for(k=0;k<N3;k++)
                    {
                        index = i*N2*N3 + j*N3 + k + plane*lsize;
                        index1 = i*N2*N3 + j*N3 + k + 0*lsize + plane*lsize*3;
                        index2 = i*N2*N3 + j*N3 + k + 1*lsize + plane*lsize*3;
                        index3 = i*N2*N3 + j*N3 + k + 2*lsize + plane*lsize*3;
                        
                        /*no derivatives taken use to calculate E_core*/
                        
                        fcore[index] += (isf/2.0)*(sin(pi*(data_core[index3].re - data_core[index1].re))*sin(pi*(data_core[index3].re - data_core[index1].re)) + sin(pi*(data_core[index1].re - data_core[index2].re))*sin(pi*(data_core[index1].re - data_core[index2].re)) + sin(pi*(data_core[index2].re - data_core[index3].re))*sin(pi*(data_core[index2].re - data_core[index3].re))) + An*(sin(2*pi*(data_core[index3].re - data_core[index1].re)+(pi*p))*sin(2*pi*(data_core[index3].re - data_core[index1].re)+(pi*p)) + sin(2*pi*(data_core[index1].re - data_core[index2].re)+(2.0*pi*p))*sin(2*pi*(data_core[index1].re - data_core[index2].re)+(2.0*pi*p)) + sin(2*pi*(data_core[index2].re - data_core[index3].re)+(pi*p))*sin(2*pi*(data_core[index2].re - data_core[index3].re)+(pi*p)))/N1;
                        
                        MPI_Reduce(&fcore[index], &f_core[index], 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
                        
                        if(rank == 0){
                            E_core += Bn*(f_core[index]/nsize);
                        }
                        df1core[index] = ((isf/2.0)*pi*(sin(2.0*pi*(data_core[index1].re-data_core[index2].re)) - sin(2.0*pi*(data_core[index3].re-data_core[index1].re)))) + An*2.0*pi*(sin(2.0*pi*(2.0*(data_core[index1].re-data_core[index2].re))+(2.0*p)) - sin(2.0*pi*(2.0*(data_core[index3].re-data_core[index1].re))+p));
                        
                        df2core[index] = ((isf/2.0)*pi*(sin(2.0*pi*(data_core[index2].re-data_core[index3].re)) - sin(2.0*pi*(data_core[index1].re-data_core[index2].re)))) + An*2.0*pi*(sin(2.0*pi*(2.0*(data_core[index2].re-data_core[index3].re)+p)) - sin(2.0*pi*(2.0*(data_core[index1].re-data_core[index2].re)+(2.0*p))));
                        
                        df3core[index] = ((isf/2.0)*pi*(sin(2.0*pi*(data_core[index3].re-data_core[index1].re)) - sin(2.0*pi*(data_core[index2].re-data_core[index3].re)))) + An*2.0*pi*(sin(2.0*pi*(2.0*(data_core[index3].re-data_core[index1].re)+p)) - sin(2.0*pi*(2.0*(data_core[index2].re-data_core[index3].re)+p)));
                    }/*ijk*/
        }/*plane*/
        
        for(plane=0;plane<NP;plane++)
        {
            for(i=0;i<lN1;i++)
                for(j=0;j<N2;j++)
                    for(k=0;k<N3;k++)
                    {
                        index = i*N2*N3 + j*N3 + k + plane*lsize;
                        index1 = i*N2*N3 + j*N3 + k + 0*lsize + plane*lsize*3;
                        index2 = i*N2*N3 + j*N3 + k + 1*lsize + plane*lsize*3;
                        index3 = i*N2*N3 + j*N3 + k + 2*lsize + plane*lsize*3;
                        
                        dE_core[index1] = Bn*df1core[index];
                        dE_core[index2] = Bn*df2core[index];
                        dE_core[index3] = Bn*df3core[index];
                    }
        }
        
        if(rank == 0 && /*fmod((double)(it),1000.0) == 0.0*/ it == NT-1){
            printf("Core Energy for this time step %lf\n", E_core);
        }
        
        MPI_Barrier(MPI_COMM_WORLD);
        //fix right now only prints first point on each processor.?????
        //if(it == NT-1){
        //if(fmod((double)(it),1000.0) == 0.0 || it == NT-1){
        if(/*fmod((double)(itp),5.0) == 0.0 ||*/ it == NT-1 && (itp == (NSI -1))){
            for(i=0;i<lsize*NP;i++){
                delta[i] = 0.0;
            }
            for(plane=0;plane<NP;plane++)
            {
                for(i=0;i<lN1;i++)
                    for(j=0;j<N2;j++)
                        for(k=0;k<N3;k++)
                        {
                            index = i*N2*N3 + j*N3 + k + plane*lsize;
                            //indexdx = ((lxs+i)+1)*N2*N3 + j*N3 + k + plane*lsize;
                            index1 = i*N2*N3 + j*N3 + k + 0*lsize + plane*lsize*3;
                            index2 = i*N2*N3 + j*N3 + k + 1*lsize + plane*lsize*3;
                            index3 = i*N2*N3 + j*N3 + k + 2*lsize + plane*lsize*3;
                            
                            delta[index] = /*(1.0/2.0)**/(data_core[index1].re*xb[0][0] + data_core[index2].re*xb[1][0] + data_core[index3].re*xb[2][0])*xb[1][0] + (data_core[index1].re*xb[0][1] + data_core[index2].re*xb[1][1] + data_core[index3].re*xb[2][1])*xb[1][1] + (data_core[index1].re*xb[0][2] + data_core[index2].re*xb[1][2] + data_core[index3].re*xb[2][2])*xb[1][2];
                            
                            //if(lxs+i == N1/2 && k == N3/2){
                            if(j == N2/2 && k == N3/2){
                                //if(rank == 0){
                                //fprintf(of4,"zone   J = %d\n", N2);
                                //}
                                fprintf(of4, "%d %lf %lf\n", lxs+i, (dx)*(lxs+i), delta[index]);
                            }
                        }/*ijk*/
            }/*plane*/
            
            for(plane=0;plane<NP;plane++)
            {
                for(i=0;i<lN1;i++)
                    //for(j=0;j<N2;j++)
                    //for(k=0;k<N3;k++)
                {
                    //if(j == N2/2 && k == N3/2){
                    j = N2/2;
                    k = N3/2;
                    index = i*N2*N3 + j*N3 + k + plane*lsize;
                    //indexdx = i*N2*N3 + (j+1)*N3 + k + plane*lsize;
                    indexdx = (i+1)*N2*N3 + j*N3 + k + plane*lsize;
                    if((i+1) == lN1){
                        indexmin = 0*N2*N3 + j*N3 + k + plane*lsize;
                        if(rank != 0){
                            MPI_Send(&delta[indexmin], 1, MPI_DOUBLE, rank-1, tag, MPI_COMM_WORLD);
                        }
                        if((lxs+(i+1)) == N1){
                            mpidel = delta[index];
                            ddelta[i] = (mpidel - delta[index])/dx;
                            goto skip_point3;
                        }
                        MPI_Recv(&mpidel, 1, MPI_DOUBLE, rank+1, tag, MPI_COMM_WORLD, &status);
                        ddelta[i] = (mpidel - delta[index])/dx;
                    }
                    else{
                        ddelta[i] = (delta[indexdx] - delta[index])/dx;
                    }
                    
                    //ddelta[index] = (delta[indexdx] - delta[index])/dx;
                    //if(lxs+i == N1/2 && k == N3/2){
                    //if(rank == 0){
                    //fprintf(of5,"zone   J = %d", N2);
                    //}
                    //fprintf(of5, "%d %lf %lf %d %d %lf %lf %lf\n", j, (dslip/N1)*j, ddelta[index], index, indexdx, delta[index], delta[indexdx], dx);
                    
                skip_point3:
                    fprintf(of5, "%d %lf %lf\n", lxs+i, (dx)*(lxs+i), ddelta[i]);
                }/*ijk*/
            }/*plane*/
        }/*it*/
    }/*cflag*/
    
    //if(fmod((double)(it),1000.0) == 0.0 || it == NT-1){
    /*if((fmod((double)(itp),5.0) == 0.0 && it == NT-1) || it == NT-1 && (itp == (NSI -1))){
     fclose(of4);
     fclose(of5);
     //fclose(of6); //--cflag == 2 only
     }*/
    /*if(rank == 0){
     fclose(of7);
     }*/
    
    free(delta);
    free(ddelta);
    free(f_core);
    
    MPI_Barrier(MPI_COMM_WORLD);
    
    return;
}

/******************************************************************/

void Bmatrix (double *BB, double *fx, double *fy, double *fz, double **xn, double **xb, double ***eps, double d1, double d2,double d3, double C[][3][3][3], double mu, double nu, double b, double dslip, int lN1, int rank, int lsize)
{
    if(rank == 0){
        printf("Set B Matrix \n");
    }
    
#define 	DELTA(i, j)   ((i==j)?1:0)
    
    int i, j, k, l, m, n, u, v, k1, k2, k3, ka, kb, nv, nb, nfreq;
    int is, js, ks;
    double fkr;
    double A[ND][ND][ND][ND];
    double B[NS][NS][lN1][N2][N3];
    double G[ND][ND];
    double fk[ND];
    double xn_temp[NS][ND], xb_temp[NS][ND], eps_temp[NS][ND][ND];
    double fk2, fk4, fka, fkb;
    
    for (i=0; i<NS; i++) {
        for (j=0; j<ND; j++) {
            xn_temp[i][j]=xn[i][j];
            xb_temp[i][j]=xb[i][j];
        }
    }
    
    /*set eps */
    
    for (i=0; i<ND; i++) {
        for (j=0; j<ND; j++) {
            for (k=0; k<NS;k++){
                if(i==j)
                    eps_temp[k][i][j] = xb_temp[k][i]*xn_temp[k][j]/dslip;
                else
                    eps_temp[k][i][j] = (xb_temp[k][i]*xn_temp[k][j] + xb_temp[k][j]*xn_temp[k][i] )/2./dslip;
            }
        }
    }
    
    for (i=0; i<ND; i++) {
        for (j=0; j<ND; j++) {
            for (k=0; k<NS;k++){
                eps[k][i][j] = eps_temp[k][i][j];
            }
        }
    }
    /* set A, Green function and B matrix*/
    for(k1=0;k1<lN1;k1++){
        for(k2=0;k2<N2;k2++){
            for(k3=0;k3<N3;k3++){
                nfreq = k3+(k2)*N3+(k1)*N3*N2;
                fk[0] = fx[nfreq];
                fk[1] = fy[nfreq];
                fk[2] = fz[nfreq];
                fk2 = fk[0]*fk[0]+fk[1]*fk[1]+fk[2]*fk[2];
                fk4 = fk2*fk2;
                
                G[0][0] = 0.0;
                
                //if(fk2>0){
                for (m=0; m<ND; m++) {
                    for (n=0; n<ND; n++) {
                        for (u=0; u<ND; u++) {
                            for (v=0; v<ND; v++) {
                                A[m][n][u][v] = 0.0;
                                
                                for	(i=0; i<ND; i++) {
                                    for (j=0; j<ND; j++) {
                                        for (k=0; k<ND; k++) {
                                            if(fk2>0)
                                                G[k][i] = (2.0 * DELTA(i,k)/fk2-1.0/(1.0-nu)*fk[i]*fk[k]/fk4)/(2.0*mu);
                                            for(l=0; l<ND; l++) {
                                                A[m][n][u][v] = A[m][n][u][v] - C[k][l][u][v]*C[i][j][m][n]*G[k][i]*fk[j]*fk[l] ;
                                            }
                                        }
                                    }
                                }
                                A[m][n][u][v] = A[m][n][u][v]+C[m][n][u][v];
                            }//v
                        }//u
                    }//n
                } //m
                
                //} /*if fk2 */
                
                
                for(ka=0;ka<NS;ka++){
                    for(kb=0;kb<NS;kb++){
                        B[ka][kb][k1][k2][k3] = 0.0;
                        for (m=0; m<ND; m++) {
                            for (n=0; n<ND; n++) {
                                for (u=0; u<ND; u++) {
                                    for (v=0; v<ND; v++) {
                                        B[ka][kb][k1][k2][k3]= B[ka][kb][k1][k2][k3] + A[m][n][u][v]*eps_temp[ka][m][n]*eps_temp[kb][u][v];
                                    }
                                }
                            }
                        }
                        
                        nb = nfreq +(ka)*lsize+(kb)*lsize*NS;
                        BB[nb] = B[ka][kb][k1][k2][k3]/mu;
                        if(k1==0&&k2==0&&k3==0){
                            BB[nb]=0.0;// Cauchy principal value
                        }
                        
                    } /*ka*/
                }/* kb*/
                
                
            }	/*k1*/
        }	/*k2*/
    }	/*k3*/
    
    return;
}

/********************************************************************/

void Fmatrix (double *FF, double *DD, double *UU, double *fx, double *fy, double *fz, double ***eps, double d1, double d2,double d3, double C[3][3][3][3], double mu, double nu, int lN1, int rank, int lsize)
{
    if(rank == 0){
        printf("set F matrix \n");
    }
    
#define 	DELTA(i, j)   ((i==j)?1:0)
    
    int i, j, k, l, m, n, u, v, k1, k2, k3, ka, nv, nb, nfreq;
    int is, js, ks;
    double fkr;
    double F[NS][ND][ND];
    double D[NS][ND][ND];
    double G[ND][ND];
    double fk[ND];
    double fk2, fk4, fka,fkb;
    double A[ND][ND][ND][ND];
    
    
    //add dis
    double U[NS][ND];
    int nc;
    /* set Green function and F matrix*/
    
    for(k1=0;k1<lN1;k1++)
        for(k2=0;k2<N2;k2++)
            for(k3=0;k3<N3;k3++)
            {
                nfreq = k3+(k2)*N3+(k1)*N3*N2;
                fk[0] = fx[nfreq];
                fk[1] = fy[nfreq];
                fk[2] = fz[nfreq];
                fk2 = fk[0]*fk[0]+fk[1]*fk[1]+fk[2]*fk[2];
                fk4 = fk2*fk2;
                
                G[0][0] = 0.0;
                for(m=0; m<ND; m++)
                    for(n=0; n<ND; n++)
                        for(u=0; u<ND; u++)
                            for(v=0; v<ND; v++)
                            {
                                A[m][n][u][v] = 0.0;
                                
                                for(i=0; i<ND; i++)
                                    for(j=0; j<ND; j++)
                                        for (k=0; k<ND; k++)
                                        {
                                            if(fk2>0)
                                                G[k][i] = (2.0 * DELTA(i,k)/fk2-1.0/(1.0-nu)*fk[i]*fk[k]/fk4)/(2.0*mu);
                                            for(l=0; l<ND; l++) 
                                            {
                                                A[m][n][u][v] = A[m][n][u][v] - C[k][l][u][v]*C[i][j][m][n]*G[k][i]*fk[j]*fk[l] ;	
                                            }	
                                        }
                                A[m][n][u][v] = A[m][n][u][v]+C[m][n][u][v];
                            }		  		
                
                for (ka=0; ka<NS; ka++)
                    for (i=0; i<ND; i++)
                        for (j=0; j<ND; j++) {
                            F[ka][i][j]= 0.0;
                            D[ka][i][j]=0.0;
                            
                            for(k=0; k<ND; k++)
                                for (l=0; l<ND; l++){
                                    for (m=0; m<ND; m++)
                                        for (n=0; n<ND; n++) {
                                            
                                            F[ka][i][j] = F[ka][i][j] + C[k][l][m][n]*G[j][k]*eps[ka][m][n]*fk[i]*fk[l] ;	
                                            
                                        }
                                    
                                    D[ka][i][j]=D[ka][i][j]+A[i][j][k][l]*eps[ka][k][l];
                                    
                                } //k,l
                            
                            nb = nfreq + (ka)*lsize + i*lsize*NS + j*lsize*NS*ND;
                            FF[nb] = F[ka][i][j];
                            DD[nb] = D[ka][i][j];
                        }  //j, i, ka	
                
                for(ka=0;ka<NS;ka++){
                    for(j=0;j<ND;j++){
                        U[ka][j] = 0.0;
                        for(k=0; k<ND; k++)
                            for (l=0; l<ND; l++)						
                                for (m=0; m<ND; m++)
                                    for (n=0; n<ND; n++){								
                                        U[ka][j] += C[k][l][m][n]*G[j][k]*eps[ka][m][n]*fk[l] ;
                                    }
                        nc = nfreq + ka*lsize + j*lsize*NS;
                        UU[nc] = U[ka][j];
                    }
                }		
                
            }/*k1,k2,k3*/
    return;
}

/********************************************************************/

void resolSS(double **sigma, double *tau, double ***eps, int rank, int ND, int NS, double C[][3][3][3], double mu)
{
    int is, i, j, k, l;
    
    for (is=0;is<NS;is++){
        tau[is] = 0.0;
        for(i=0;i<ND;i++){
            for(j=0;j<ND;j++){
                if(load == 0)			//stress controlled
                    tau[is] = tau[is]+(sigma[i][j]*eps[is][i][j]);
                else {					//strain controlled
                    for(k=0;k<ND;k++)
                        for(l=0;l<ND;l++){
                            tau[is] += C[i][j][k][l]*epsilon[i][j]*eps[is][k][l]/mu;
                        }
                }
            }	
        }
        if(rank == 0){
            printf("Resolved shear stresses tau[%d] = %lf\n", is, tau[is]);
        }
    }
    
    return;
}

/********************************************************************/

void material(int mflag, double *a, double *mu, double *young, double *c0, double *c1, double *c2, double *c3, double *c4, double *a1, double *a3, double *isf, double *usf, double *nu, double *ll, double *C44, double *C12, double *C11, double *S44, double *S12, double *S11, double C[][3][3][3])
{
    /*mflag == 1 -> Nickel (PRISM device)
     mflag == 2 -> Aluminum
     mflag == 3 -> Gold
     mflag == 4 -> Copper
     mflag == 5 -> Acetaminophen
     mflag == 6 -> Palladium*/
    
    if(mflag == 1){ /*Nickel*/
        *a = 3.52E-10; //meters
        *mu = 75.0E9; //Pascal = N/m^2 = kg/(ms^2)
        *young = 200.0E9; //Pascal = N/m^2 = kg/(ms^2)
        
        *c0 = 410.024E-3; //J/m^2
        *c1 = -51.99716E-3; //gamma surface coefficients
        *c2 = -120.555E-3;
        *c3 = 35.21191E-3;
        *c4 = 0.594646E-3;
        *a1 = -66.1927E-3;
        *a3 = -75.3124E-3;
        //coefficients from Strachan etal MD simulations
        
        *isf = 84.718E-3; //intrinsic stacking fault energy
        *usf = 211.688E-3; //unstable stacking fault energy
        //J/m^2 
    }
    if(mflag == 2){ /*Aluminum*/
        *a = 4.05E-10; //meters
        *mu = 26.0E9; //Pascal
        *young = 70.0E9; //Pascal
        
        *c0 = 242.5E-3; //J/m^2
        *c1 = -51.65E-3; //gamma surface coefficients
        *c2 = -39.71E-3;
        *c3 = 13.68E-3;
        *c4 = -1.572E-3;
        *a1 = -30.40E-3;
        *a3 = -13.72E-3;
        //gamma-surface coefficients from Shen and Wang, 2004 units mJ/m2
        
        *isf = 141.786E-3;  //intrinsic stacking fault energy 
        *usf = 385.3457E-3; //unstable stacking fault energy 
        //J/m2 simulated gamma surface via Shen and Wang 2004
    }
    if(mflag == 3){ /*Gold*/
        *a = 4.08E-10; //meters
        *mu = 27.0E9; //Pascal 
        *young = 78.0E9; //Pascal 
        
        *c0 = 0.0; //J/m^2
        *c1 = 0.0; //gamma surface coefficients
        *c2 = 0.0;
        *c3 = 0.0;
        *c4 = 0.0;
        *a1 = 0.0;
        *a3 = 0.0;
        
        *isf = 37.0E-3; //intrinsic stacking fault energy
        *usf = 110.0E-3; //unstable stacking fault energy
        //J/m2 from Tadmor and Bernstein 2004 
    }
    if(mflag == 4){ /*Copper*/
        *a = 3.61E-10; //meters
        *mu = 48.0E9; //Pascal 
        *young = 110.0E9; //Pascal
        
        *c0 = 0.0; //J/m^2
        *c1 = 0.0; //gamma surface coefficients
        *c2 = 0.0;
        *c3 = 0.0;
        *c4 = 0.0;
        *a1 = 0.0;
        *a3 = 0.0;
        
        *isf = 0.0; //intrinsic stacking fault energy
        *usf = 0.0; //unstable stacking fault energy
    }
    if(mflag == 5){ /*Acetaminophen*/
        *a = 1.84E-9; //meters
        *mu = 3.2E9; //Pascal
        *young = 8.0E9; //Pascal
        
        *c0 = 0.0; //J/m^2
        *c1 = 0.0; //gamma surface coefficients
        *c2 = 0.0;
        *c3 = 0.0;
        *c4 = 0.0;
        *a1 = 0.0;
        *a3 = 0.0;
        
        *isf = 0.0; //intrinsic stacking fault energy
        *usf = 0.0; //unstable stacking fault energy
    }
    if(mflag == 6){ /*Palladium*/
        *a = 3.890E-10; //meters
        *mu = 44.0E9;  //Pascal 
        *young = 121.0E9; //Pascal
        
        *c0 = 374.1E-3; //J/m^2
        *c1 = -71.01E-3; //gamma surface coefficients
        *c2 = -69.62E-3;
        *c3 = 21.22E-3;
        *c4 = -2.644E-3;
        *a1 = -54.80E-3;
        *a3 = -27.13E-3;
        //gamma-surface coefficients from Shen and Wang, 2004 units mJ/m2
        
        *isf = 177.8237E-3; //intrinsic stacking fault energy 
        *usf = 596.9907E-3; //unstable stacking fault energy
        //J/m2 simulated gamma surface via Shen and Wang 2004
    }
    
    *C44 = *mu; //Pa
    *nu = *young/2.0/(*mu)-1.0;
    *C12 = 2.0**nu**C44/(1.0-2.0**nu); //Pa
    *C11 = 2.0**C44+*C12;//Pa
    *ll = *C12;
    
    *S11 = 1.0/(*young);
    *S12 = -*nu/(*young);
    *S44 = 2*(*S11-*S12);
    
    int i, j, k, m;
    double mup= *C11-*C12-2*(*C44);
    
#define 	DELTA(i, j)   ((i==j)?1:0)
#define		DELTA4(i,j,k,l) (((i==j) && (j==k) && (k==l))?1:0)
    
    for (i=0; i<ND; i++) 
        for (j=0; j<ND; j++) 
            for (k=0; k<ND; k++)
                for (m=0; m<ND; m++) {
                    C[i][j][k][m] = (*mu)*(DELTA(i,k)*DELTA(j,m)+DELTA(i,m)*DELTA(j,k)) + (*ll)*DELTA(i,j)*DELTA(k,m) + mup*DELTA4(i,j,k,m);	
                }
    return;
}

/********************************************************************/

void set2D(double **xn, double **xb, int rank, int oflag)
{
    xn[0][0]= 0.0; //1.0/sqrt(3);
    xn[0][1]= 0.0; //1.0/sqrt(3);
    xn[0][2]= 1.0; //1.0/sqrt(3);
    
    if(oflag == 0){
        xb[0][0]= 1.0; //-1.0/2.0;
        xb[0][1]= 0.0; //1.0/2.0;
        xb[0][2]= 0.0;
    }
    else if(oflag == 1){
        xb[0][0]= 0.0; //-1.0/2.0;
        xb[0][1]= 1.0; //1.0/2.0;
        xb[0][2]= 0.0;
    }
    else{
        printf("Initial configuration specified has not been developed yet"); 
    }
    
    if(rank == 0){
        printf("Burgers vector b = (%lf , %lf , %lf )\n", xb[0][0],xb[0][1], xb[0][2]);	
        printf("Slip plane     n = [%lf,  %lf , %lf ]\n", xn[0][0],xn[0][1], xn[0][2]);
    }
    
    return;
}

/********************************************************************/

void set3D1pl(double **xn, double **xb, int rank, int oflag)
{
    
    int i, j, k;
    double s3;
    
    s3 = sqrt(3);
    for (i=0;i<3; i++) {
        xn[i][0]= 0.0; //1.0/s3;
        xn[i][1]= 0.0; //1.0/s3;
        xn[i][2]= 1.0; //1.0/s3;
    }
    
    if(oflag == 0){
        xb[0][0]= -0.5; //0.0; //1.0;
        xb[0][1]= -sqrt(3.0)/2.0; //1.0; //0.0;
        xb[0][2]= 0.0; //0.0; //0.0;
        
        xb[1][0]= 1.0; //sqrt(3.0)/2.0; //-0.5 
        xb[1][1]= 0.0; //-0.5; //sqrt(3.0)/2.0;
        xb[1][2]= 0.0; //0.0; //0.0;
        
        xb[2][0]= -0.5;//-sqrt(3.0)/2.0; //-0.5;
        xb[2][1]= sqrt(3.0)/2.0; //-0.5; //-sqrt(3.0)/2.0;
        xb[2][2]= 0.0; //0.0; //0.0;
    }
    else if(oflag == 1){
        xb[0][0]= sqrt(3.0)/2.0; //0.0; //1.0;
        xb[0][1]= -0.5; //1.0; //0.0;
        xb[0][2]= 0.0; //0.0; //0.0;
        
        xb[1][0]= 0.0; //sqrt(3.0)/2.0; //-0.5 
        xb[1][1]= 1.0; //-0.5; //sqrt(3.0)/2.0;
        xb[1][2]= 0.0; //0.0; //0.0;
        
        xb[2][0]= -sqrt(3.0)/2.0;//-sqrt(3.0)/2.0; //-0.5;
        xb[2][1]= -0.5; //-0.5; //-sqrt(3.0)/2.0;
        xb[2][2]= 0.0; //0.0; //0.0;
    }
    else{
        printf("Initial configuration specified has not been developed yet"); 
    }
    
    if(rank == 0){
        printf("Burgers vector b1 = (%lf , %lf , %lf )\n", xb[0][0],xb[0][1], xb[0][2]);
        printf("Burgers vector b2 = (%lf , %lf , %lf )\n", xb[1][0],xb[1][1], xb[1][2]);
        printf("Burgers vector b3 = (%lf , %lf , %lf )\n", xb[2][0],xb[2][1], xb[2][2]);	
        printf("Slip plane     n = [%lf,  %lf , %lf ]\n", xn[0][0],xn[0][1], xn[0][2]);
    }
    
    return;
}

/***********************************************************************/

void set3D2sys(double **xn, double **xb, int rank)
{
    xn[0][0]=1.0/sqrt(3);
    xn[0][1]=1.0/sqrt(3);
    xn[0][2]=1.0/sqrt(3);
    
    xb[0][0]=-1.0/2.0;
    xb[0][1]=1.0/2.0;
    xb[0][2]=0.0;
    
    xn[1][0]=1.0/sqrt(3);
    xn[1][1]=1.0/sqrt(3);
    xn[1][2]=1.0/sqrt(3);
    
    xb[1][0]=0.0;
    xb[1][1]=-1.0/2.0;
    xb[1][2]=1.0/2.0;
    
    if(rank == 0){
        printf("Burgers vector b1 = (%lf , %lf , %lf )\n", xb[0][0],xb[0][1], xb[0][2]);
        printf("Burgers vector b2 = (%lf , %lf , %lf )\n", xb[1][0],xb[1][1], xb[1][2]);
        printf("Slip plane     n  = [%lf,  %lf , %lf ]\n", xn[0][0],xn[0][1], xn[0][2]);
    }
    
    return;
}

/**********************************************************************/

void setfcc (double **xn, double **xb, int rank)
{
    int i, j, k;
    double s3;
    
    s3 = sqrt(3);
    for (i=0;i<3 ; i++) {
        xn[i][0]= 1.0/s3;
        xn[i][1]= 1.0/s3;
        xn[i][2]= 1.0/s3;
    }
    
    for (i=3;i<6 ; i++) {
        xn[i][0]= -1.0/s3;
        xn[i][1]= 1.0/s3;
        xn[i][2]= 1.0/s3;
    }
    
    for (i=6;i<9 ; i++) {
        xn[i][0]= 1.0/s3;
        xn[i][1]= -1.0/s3;
        xn[i][2]= 1.0/s3;
    }
    for (i=9;i<12 ; i++) {
        xn[i][0]= 1.0/s3;
        xn[i][1]= 1.0/s3;
        xn[i][2]= -1.0/s3;
    }
    xb[0][0] = -1.0;
    xb[0][1] = 1.0;
    xb[0][2] = 0.0;
    
    xb[1][0] = 0.0;
    xb[1][1] = -1.0;
    xb[1][2] = 1.0;
    
    xb[2][0] = 1.0;
    xb[2][1] = 0.0;
    xb[2][2] = -1.0;
    
    xb[3][0] = -1.0;
    xb[3][1] = -1.0;
    xb[3][2] = 0.0;
    
    xb[4][0] = 1.0;
    xb[4][1] = 0.0;
    xb[4][2] = 1.0;
    
    xb[5][0] = 0.0;
    xb[5][1] = -1.0;
    xb[5][2] = 1.0;
    
    xb[6][0] = -1.0;
    xb[6][1] = -1.0;
    xb[6][2] = 0.0;
    
    xb[7][0] = 1.0;
    xb[7][1] = 0.0;
    xb[7][2] = -1.0;
    
    xb[8][0] = 0.0;
    xb[8][1] = -1.0;
    xb[8][2] = -1.0;
    
    xb[9][0] = -1.0;
    xb[9][1] = 1.0;
    xb[9][2] = 0.0;
    
    xb[10][0] = 1.0;
    xb[10][1] = 0.0;
    xb[10][2] = 1.0;
    
    xb[11][0] = 0.0;
    xb[11][1] = -1.0;
    xb[11][2] = -1.0;
    
    if(rank == 0){
        printf("Slip systems\n");
    }
    for (i=0; i<12; i++){
        for (j=0; j<3; j++) {
            xb[i][j] = xb[i][j]/2.0;	
        }
        if(rank == 0){
            printf("b(%d) = (%lf , %lf , %lf )\n", i, xb[i][0],xb[i][1], xb[i][2]);	
            printf("n(%d) = [%lf,  %lf , %lf ]\n", i, xn[i][0],xn[i][1], xn[i][2]);
        }   
    }
    return;	
}

/********************************************************************/

void frec( double *fx,double *fy, double *fz, double d1, double d2, double d3,
          int lN1, int lxs, int N1, int N2, int N3)
{
    int i,j,k,ksym, nf; 
    
    for(i=0;i<lN1;i++)  //need to pass ln1
    {
        for(j=0;j<N2;j++)
        {
            for(k=0;k<N3;k++)
            {
                nf = k+(j)*N3+(i)*N3*N2;
                /* frecuency in x */
                if (lxs+i==0) {  //need to pass lxs
                    fx[nf]= 0.0;
                }
                if (lxs+i >= 1 && lxs+i < N1/2 ) {
                    fx[nf]= (double)(lxs+i)/((double)(N1)/d1);
                }
                if (lxs+i >= N1/2) {
                    fx[nf]= ((double)(lxs+i)-(double)(N1))/(double)(N1)/d1;	
                }
                /* frecuency in y */
                if (j==0) {
                    fy[nf]= 0.0;
                }
                if (j >= 1 && j < N2/2 ) {
                    fy[nf]= (double)(j)/(double)(N2)/d2;
                }
                if (j >= N2/2) {
                    fy[nf]= ((double)(j)-(double)(N2))/(double)(N2)/d2;
                }				
                /* frecuency in z */
                if (k==0) {
                    fz[nf]= 0.0;
                }
                if (k >= 1 && k < N3/2 ) {
                    fz[nf]= (double)(k)/(double)(N3)/d3;
                }	
                if (k >= N3/2) {
                    fz[nf]= ((double)(k)-(double)(N3))/(double)(N3)/d3;
                }		
                
                /*printf("%d %d %d    %lf %lf %lf \n", i, j, k, fx[nf], fy[nf],fz[nf]); 	*/
                
            }
        }
    } 
    return;
}

/********************************************************************/

/*void seteps (double eps[NS][ND][ND], double epsv[NV][ND][ND],double xn[NS][ND], double xb[NS][ND], double dslip)
 {
 double xn1[NV][ND], xb1[NV][ND];
 int i, j, k, v;
	
 for (v=0;v<NV;v++){
 for(i=0;i<ND;i++){
 xn1[v][i]=0.0;
 xb1[v][i]=0.0;
 }	
 }
	
 for (v=0;v<NV;v++){
 if(v==0 || v==1 || v==2)
 {xn1[v][0]=1.0;}
 if(v==3 || v==4 || v==5)
 {xn1[v][1]=1.0;}
 if(v==6 || v==7 || v==8)
 {xn1[v][2]=1.0;}
 if(v==0 || v==3 || v==6)
 {xb1[v][0]=1.0;}
 if(v==1 || v==4 || v==7)
 {xb1[v][1]=1.0;}
 if(v==2 || v==5 || v==8)
 {xb1[v][2]=1.0;}
 // printf("%lf %lf %lf %lf %lf %lf \n", xn1[v][0],xn1[v][1],xn1[v][2],xb1[v][0],xb1[v][1],xb1[v][2]);
 }
 
 for (v=0; v<NV;v++){
 printf("strainsys (%d) \n", v);	
 for (i=0; i<ND; i++) {
 for (j=0; j<ND; j++) {
	epsv[v][i][j]= xn1[v][i]*xb1[v][j];	
	printf("%lf ", epsv[v][i][j]);				
 }
 printf("\n");
 }
 }
 
 /*set eps */

/*for (i=0; i<ND; i++) {
 for (j=0; j<ND; j++) {
 for (k=0; k<NS;k++){
	//eps[k][i][j]= xb[k][i]*xn[k][j]/dslip;
	if(i==j)
	eps[k][i][j] = xb[k][i]*xn[k][j]/dslip;
	else
 eps[k][i][j] = (xb[k][i]*xn[k][j] + xb[k][j]*xn[k][i] )/dslip;				
 }
 }
 }
	
 return;
 }*/ 

/********************************************************************/

/* 
 A C-program for MT19937, with initialization improved 2002/1/26.
 Coded by Takuji Nishimura and Makoto Matsumoto.
 
 Before using, initialize the state by using init_genrand(seed)  
 or init_by_array(init_key, key_length).
 
 Copyright (C) 1997 - 2002, Makoto Matsumoto and Takuji Nishimura,
 All rights reserved.                          
 
 Redistribution and use in source and binary forms, with or without
 modification, are permitted provided that the following conditions
 are met:
 
 1. Redistributions of source code must retain the above copyright
 notice, this list of conditions and the following disclaimer.
 
 2. Redistributions in binary form must reproduce the above copyright
 notice, this list of conditions and the following disclaimer in the
 documentation and/or other materials provided with the distribution.
 
 3. The names of its contributors may not be used to endorse or promote 
 products derived from this software without specific prior written 
 permission.
 
 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR
 CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 
 
 Any feedback is very welcome.
 http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/emt.html
 email: m-mat @ math.sci.hiroshima-u.ac.jp (remove space)
 */

#include <stdio.h>

/* Period parameters */  
#define N 624
#define M 397
#define MATRIX_A 0x9908b0dfUL   /* constant vector a */
#define UPPER_MASK 0x80000000UL /* most significant w-r bits */
#define LOWER_MASK 0x7fffffffUL /* least significant r bits */

static unsigned long mt[N]; /* the array for the state vector  */
static int mti=N+1; /* mti==N+1 means mt[N] is not initialized */

/* initializes mt[N] with a seed */
void init_genrand(unsigned long s)
{
    mt[0]= s & 0xffffffffUL;
    for (mti=1; mti<N; mti++) {
        mt[mti] = 
        (1812433253UL * (mt[mti-1] ^ (mt[mti-1] >> 30)) + mti); 
        /* See Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier. */
        /* In the previous versions, MSBs of the seed affect   */
        /* only MSBs of the array mt[].                        */
        /* 2002/01/09 modified by Makoto Matsumoto             */
        mt[mti] &= 0xffffffffUL;
        /* for >32 bit machines */
    }
}

/* initialize by an array with array-length */
/* init_key is the array for initializing keys */
/* key_length is its length */
/* slight change for C++, 2004/2/26 */
void init_by_array(unsigned long init_key[], int key_length)
{
    int i, j, k;
    init_genrand(19650218UL);
    i=1; j=0;
    k = (N>key_length ? N : key_length);
    for (; k; k--) {
        mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 30)) * 1664525UL))
        + init_key[j] + j; /* non linear */
        mt[i] &= 0xffffffffUL; /* for WORDSIZE > 32 machines */
        i++; j++;
        if (i>=N) { mt[0] = mt[N-1]; i=1; }
        if (j>=key_length) j=0;
    }
    for (k=N-1; k; k--) {
        mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 30)) * 1566083941UL))
        - i; /* non linear */
        mt[i] &= 0xffffffffUL; /* for WORDSIZE > 32 machines */
        i++;
        if (i>=N) { mt[0] = mt[N-1]; i=1; }
    }
    
    mt[0] = 0x80000000UL; /* MSB is 1; assuring non-zero initial array */ 
}

/* generates a random number on [0,0xffffffff]-interval */
unsigned long genrand_int32(void)
{
    unsigned long y;
    static unsigned long mag01[2]={0x0UL, MATRIX_A};
    /* mag01[x] = x * MATRIX_A  for x=0,1 */
    
    if (mti >= N) { /* generate N words at one time */
        int kk;
        
        if (mti == N+1)   /* if init_genrand() has not been called, */
            init_genrand(5489UL); /* a default initial seed is used */
        
        for (kk=0;kk<N-M;kk++) {
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk+M] ^ (y >> 1) ^ mag01[y & 0x1UL];
        }
        for (;kk<N-1;kk++) {
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk+(M-N)] ^ (y >> 1) ^ mag01[y & 0x1UL];
        }
        y = (mt[N-1]&UPPER_MASK)|(mt[0]&LOWER_MASK);
        mt[N-1] = mt[M-1] ^ (y >> 1) ^ mag01[y & 0x1UL];
        
        mti = 0;
    }
    
    y = mt[mti++];
    
    /* Tempering */
    y ^= (y >> 11);
    y ^= (y << 7) & 0x9d2c5680UL;
    y ^= (y << 15) & 0xefc60000UL;
    y ^= (y >> 18);
    
    return y;
}

/* generates a random number on [0,0x7fffffff]-interval */
long genrand_int31(void)
{
    return (long)(genrand_int32()>>1);
}

/* generates a random number on [0,1]-real-interval */
double genrand_real1(void)
{
    return genrand_int32()*(1.0/4294967295.0); 
    /* divided by 2^32-1 */ 
}

/* generates a random number on [0,1)-real-interval */
double genrand_real2(void)
{
    return genrand_int32()*(1.0/4294967296.0); 
    /* divided by 2^32 */
}

/* generates a random number on (0,1)-real-interval */
double genrand_real3(void)
{
    return (((double)genrand_int32()) + 0.5)*(1.0/4294967296.0); 
    /* divided by 2^32 */
}

/* generates a random number on [0,1) with 53-bit resolution*/
double genrand_res53(void) 
{ 
    unsigned long a=genrand_int32()>>5, b=genrand_int32()>>6; 
    return(a*67108864.0+b)*(1.0/9007199254740992.0); 
} 

void inputVariables(char *PFDDinput){
    // Open the file for reading
    FILE *fp;
    fp = fopen(PFDDinput,"r");
    if(fp == NULL){
        printf("ERROR READING THE FILE, MAYBE THE FILE DON'T EXIST");
    }
    char line[200];
    int linenum=0;
    char value[30];
    char variable [50]; 	
    
    // Reads the values from the input file(fp) to the variable 'line' line by line
    
    while(fgets(line, 200, fp) != NULL){
        linenum++;
        /** If the first character in the line is an '#' is because it's a comment 
         * so it will not take that line, it will continue with the next line
         */
        if(line[0] == '#' || line[0]== ' '){ 
            continue;
        }
        value[0] = '\0';
        variable[0] ='\0';
        // This part is the one who gives truble, hope i have fix it
        if(sscanf(line, "%s %s",variable, value) != 2){ 
            // Check if get to the end of the file (IT IS NEEDED????)
            if(variable[0] == '\0'){
                //This printf can be erase or commented, is used to check when we reach the end of the file
                //printf("END OF THE INPUT FILE\n");
                continue;
            } 
            fprintf(stderr, "Syntax error, line %d in inputfile where variable %s is equal to %s\n", linenum,variable,value);
            continue;
        }
        
        //For each variable set the new values it only works with if-else
        if(strcmp(variable,"N1") == 0)
            N1 = atoi(value);
        else if(strcmp(variable,"N2") == 0)
            N2 = atoi(value);			
        else if(strcmp(variable,"N3") == 0)
            N3 = atoi(value);			
        else if(strcmp(variable,"NT") == 0)
            NT = atoi(value);
        else if(strcmp(variable,"load") == 0)
            load = atoi(value);			
        else if(strcmp(variable,"NSI") == 0)
            NSI = atoi(value);			
        else if(strcmp(variable,"itp") == 0)
            itp = atoi(value);			
        else if(strcmp(variable,"rsm") == 0)
            rsm = atoi(value);			
        else if(strcmp(variable,"NS") == 0)
            NS = atoi(value);			
        else if(strcmp(variable,"ND") == 0)
            ND = atoi(value);			
        else if(strcmp(variable,"NP") == 0)
            NP = atoi(value);			
        else if(strcmp(variable,"iflag") == 0)
            iflag = atoi(value);			
        else if(strcmp(variable,"nloops") == 0)
            nloops = atoi(value);			
        else if(strcmp(variable,"config") == 0)
            config = atoi(value);			
        else if(strcmp(variable,"radius2") == 0)
            radius2 = atoi(value);			
        else if(strcmp(variable,"mflag") == 0)
            mflag = atoi(value);			
        else if(strcmp(variable,"cflag") == 0)
            cflag = atoi(value);			
        else if(strcmp(variable,"An") == 0)
            An = atof(value);			
        else if(strcmp(variable,"Bn") == 0)
            Bn= atof(value);			
        else if(strcmp(variable,"oflag") == 0)
            oflag = atoi(value);			
        else if(strcmp(variable,"setsigma") == 0)
            setsigma = atof(value);			
        else if(strcmp(variable,"sigstep") == 0)
            sigstep = atof(value);			
        else if(strcmp(variable,"seteps") == 0)
            seteps = atof(value);			
        else if(strcmp(variable,"epstep") == 0)
            epstep = atof(value);				
        else if(strcmp(variable,"scratch") == 0)
            strcpy(scratch, value);			
        else if(strcmp(variable,"voro") == 0)
            voro = atoi(value);			
        else if(strcmp(variable,"test") == 0)
            test = atoi(value);
        else if(strcmp(variable,"dump_3") == 0)
            dump_3 = atoi(value);
        else if(strcmp(variable,"comment") == 0)
            comment = atoi(value);
        
        // Display that the variable enter is not recongnized by the program
        else
            printf("The variable %s in the inputfile is not recognized\n",variable);
        //	printf"Line %d:  Variable  %s con el valor %d\n", linenum, variable, value);
        //	linenum++;
    }
    fclose(fp);
    
    if(test == 1)
        strcpy(StandardFile, "Standard_ss");
    
    checkVariables();
}

void checkVariables(){
    char wrong = 0;
    if (N1 == -200){
        printf("The variable N1 is not defined in the inputfile, please define\n");
        wrong = 1;
    }
    if (N2 == -200){
        printf("The variable N2 is not defined in the inputfile, please define\n");
        wrong = 1;
    }
    if (N3 == -200){
        printf("The variable N3 is not defined in the inputfile, please define\n");
        wrong = 1;
    }
    if (NT == -200){
        printf("The variable NT is not defined in the inputfile, please define\n");
        wrong = 1;
    }
    if (load == -200){
        printf("The variable load is not defined in the inputfile, please define\n");
        wrong = 1;
    }
    if (NSI == -200){
        printf("The variable NSI is not defined in the inputfile, please define\n");
        wrong = 1;
    }
    if (itp == -200){
        printf("The variable itp is not defined in the inputfile, please define\n");
        wrong = 1;
    }
    if (rsm == -200){
        printf("The variable rsm is not defined in the inputfile, please define\n");
        wrong = 1;
    }
    if (NS == -200){
        printf("The variable NS is not defined in the inputfile, please define\n");
        wrong = 1;
    }
    if (ND == -200){
        printf("The variable ND is not defined in the inputfile, please define\n");
        wrong = 1;
    }
    if (NP == -200){
        printf("The variable NP is not defined in the inputfile, please define\n");
        wrong = 1;
    }
    if (iflag == -200){
        printf("The variable iflag is not defined in the inputfile, please define\n");
        wrong = 1;
    }
    if (nloops == -200){
        printf("The variable nloops is not defined in the inputfile, please define\n");
        wrong = 1;
    }
    if (config == -200){
        printf("The variable config is not defined in the inputfile, please define\n");
        wrong = 1;
    }
    if (radius2 == -200){
        printf("The variable radius2 is not defined in the inputfile, please define\n");
        wrong = 1;
    }
    if (mflag == -200){
        printf("The variable mflag is not defined in the inputfile, please define\n");
        wrong = 1;
    }
    if (cflag == -200){
        printf("The variable cflag is not defined in the inputfile, please define\n");
        wrong = 1;
    }
    if (An == -200){
        printf("The variable An is not defined in the inputfile, please define\n");
        wrong = 1;
    }
    if (Bn == -200){
        printf("The variable Bn is not defined in the inputfile, please define\n");
        wrong = 1;
    }	
    if (oflag == -200){
        printf("The variable oflag is not defined in the inputfile, please define\n");
        wrong = 1;
    }	
    if (setsigma == -200){
        printf("The variable setsigma is not defined in the inputfile, please define\n");
        wrong = 1;
    }
    if (sigstep == -200){
        printf("The variable sigstep is not defined in the inputfile, please define\n");
        wrong = 1;
    }	
    if (seteps == -200){
        printf("The variable setsigma is not defined in the inputfile, please define\n");
        wrong = 1;
    }	
    if (epstep == -200){
        printf("The variable sigstep is not defined in the inputfile, please define\n");
        wrong = 1;
    }		
    if (voro == -200){
        printf("The variable voro is not defined in the inputfile, please define\n");
        wrong = 1;
    }
    if (scratch[0] == '\0'){
        printf("The variable scratch is not defined in the inputfile, please define\n");
        wrong = 1;
    }	
    if (test == -200){
        printf("The variable test is not defined in the inputfile, please define\n");
        wrong = 1;
    }
    if (dump_3 == -200){
        printf("The variable dump_3 is not defined in the inputfile, please define\n");
        wrong = 1;
    }
    if (comment == -200){
        printf("The variable comment is not defined in the inputfile, please define\n");
        wrong = 1;
    }
    if(wrong){
        exit(1);
    }
}

void verification(int test, double E_elas, double E_elas_test){
    if(test == 1){
        FILE *infStandard, *DiffFile;
        char line[15];
        line [0]= 'C';
        if(test == 1)
            infStandard = fopen(StandardFile,"r");
        if(infStandard == NULL){
            printf("ERROR IN THE SAMPLE FILE, DON'T EXIST OR IS NOT NAMED %s \n", StandardFile);
        }else{
            system("cat outsscurve* > sscurve.dat");
            system("diff sscurve.dat Standard_ss > Verification_ss.txt");
            fclose(infStandard);
        }
        DiffFile = fopen("Verification_ss.txt","r");
        fgets(line, 15, DiffFile);
        if(line[0] == 'C'){
            printf("\n----------------\n The verification is successful. stress-strain curve is correct.\n");
        }else{
            printf("ERROR: \n----------------\n stress-strain curve wrong. Please check Verification_ss.txt for more details.\n");
        }
        
        fclose(DiffFile);
    } // test1
    
    else if(test == 2){
        if( fabs(E_elas-3.817089e-12) <0.01e-12 && fabs((E_elas-E_elas_test)/E_elas)<1.0E-6 )
            printf("\n----------------\n The verification is successful. elastic energy is correct.\n");
        else
            printf("ERROR: \n----------------\n E_elas=%e E_elas_test=%e Standard E_elas=3.817089e-12", E_elas, E_elas_test);
    }//test2
}




/* These real versions are due to Isaku Wada, 2002/01/09 added */

/*int main(void)
 {
 int i;
 unsigned long init[4]={0x123, 0x234, 0x345, 0x456}, length=4;
 init_by_array(init, length);
 printf("1000 outputs of genrand_int32()\n");
 for (i=0; i<1000; i++) {
 printf("%10lu ", genrand_int32());
 if (i%5==4) printf("\n");
 }
 printf("\n1000 outputs of genrand_real2()\n");
 for (i=0; i<1000; i++) {
 printf("%10.8f ", genrand_real2());
 if (i%5==4) printf("\n");
 }
 return 0;
 }*/
