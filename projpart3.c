//========================================================================================================
// USING METROPOLIS-HASTINGS TO FIND MEAN SEPARATION OF ELECTRONS AND VARIANCE OF ENERGY AT ENERGY MINIMUM
//========================================================================================================
 
//note: manually change w between 0.01, 0.5 and 1.0 and run the programme
#include <stdio.h>
#include <math.h>
#include <stdint.h>
#include <time.h>
#include "vmc.h"
 
#define N_al 200 //number of steps for alpha
#define N_r 300000 //number of steps for the coordinates
#define N_f 5000000 //number of steps for final evaluation
#define B (N_r/5) //burn-in
#define PLT 25000 //number of loops for stepfinder
 
static double wval[3]={0.01, 0.5, 1.0}; double a[N_al]={0.0}, Ealpha[N_al]={0.0}, VarEalpha[N_al]={0.0}, vals[3]={0.0}; long long acpt; //vals stores the results of E, Var(E), <r12>
//acpt counts the accepted moves after burn-in
 
void eval(double al, double w, long long N){ //calculates E, Var(E), <r12> and adds the accepted moves (after burn-in) to acpt
    double step=stepfinderT1(al, w, PLT), En=0.0, En2=0.0, mr=0.0; long long bn=N/5;
    MH k; MH_init(&k);
    for(long long j=0; j<N; j++){
        int mv=MH_stepT1(&k, step, al, w); if(j<bn) continue;
        acpt+=mv;
        double e=EL1(k.r1, k.r2, k.r12, al, w);
        En+=e; En2+=e*e; mr+=k.r12;
    }
    vals[0]=En/(N-bn); vals[1]=En2/(N-bn)-vals[0]*vals[0]; vals[2]=mr/(N-bn);
}
 
int main(){
    rng_seed((uint64_t)time(NULL) ^ ((uint64_t)clock() << 32)); //initialising the RNG
    
    FILE *fp=fopen("data.csv", "w"); if(!fp) return 1;

    for(int j=0; j<3; j++){
        double w=wval[j]; acpt=0;
        double stepa=0.7/(N_al-1); //taking a uniformly over [0.3,1.0]
        for(int i=0; i<N_al; i++){
            a[i]=0.3+i*stepa; //incrementing to the next alpha with fixed interval
            eval(a[i], w, N_r);
            Ealpha[i]=vals[0]; VarEalpha[i]=vals[1]; //storing the values of E(alpha) and its variance
        }
        int idx=0;
        for(int i=0; i<N_al; i++){
            if(Ealpha[i]<Ealpha[idx]) idx=i; //finding minimum energy
            fprintf(fp, "%f,%.8e,%.8e,%.2f\n", a[i], Ealpha[i], VarEalpha[i], w);
        }
        printf("Acceptance percentage is %f\n", 100.0*acpt/((double)N_al*(N_r-B)));
        printf("Minimum E(alpha) is: %f\nCorresponding to alpha: %f\n", Ealpha[idx], a[idx]);
    
        //Computing expectation value of mean distance between electrons (and variance of energy) at energy minimum:
        double al=a[idx];
        eval(al, w, N_f); //longer run at the minimum
        printf("Expectation of mean separation for\nw=%.2f and alpha=%f is:\n<r12> = %f\nVar(E) at the minimum is: %f\n\n", w, al, vals[2], vals[1]);
    }
    fclose(fp);
    
    return 0;
}