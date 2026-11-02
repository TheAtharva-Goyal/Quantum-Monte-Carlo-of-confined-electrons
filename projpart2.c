//=================================================================================
// USING METROPOLIS-HASTINGS ALGORITHM TO PLOT ENERGY VS ALPHA FOR TRIAL FUNCTION 1
//=================================================================================
 
#include <stdio.h>
#include <math.h>
#include <stdint.h>
#include <time.h>
#include "vmc.h"
 
#define N_al 200 //number of steps for alpha
#define N_r 200000 //number of steps for the coordinates
#define B (N_r/5) //burn-in
#define PLT 25000 //number of loops for stepfinder
 
static double wval[3]={0.01, 0.5, 1.0}; double a[N_al]={0.0}, Ealpha[N_al]={0.0}; long long acpt; //acpt counts the accepted moves after burn-in
 
double eval(double al, double w){ //calculates E(alpha) and adds the accepted moves (after burn-in) to acpt
    double step=stepfinderT1(al, w, PLT), En=0.0;
    MH k; MH_init(&k);
    for(int j=0; j<N_r; j++){
        int mv=MH_stepT1(&k, step, al, w); if(j<B) continue;
        acpt+=mv;
        En+=EL1(k.r1, k.r2, k.r12, al, w);
    }
    return En/(N_r-B);
}
 
int main(){
    rng_seed((uint64_t)time(NULL) ^ ((uint64_t)clock() << 32)); //initialising the RNG
    
    FILE *fp=fopen("data.csv", "w"); if(!fp) return 1;

    for(int j=0; j<3; j++){
        double w=wval[j]; acpt=0;
        double stepa=0.9/(N_al-1); //taking a uniformly over [0.3,1.2]
        for(int i=0; i<N_al; i++){
            a[i]=0.3+i*stepa; //incrementing to the next alpha with fixed interval
            Ealpha[i]=eval(a[i], w); //storing the value of E(alpha)
        }
        int idx=0;
        for(int i=0; i<N_al; i++){
            if(Ealpha[i]<Ealpha[idx]) idx=i; //finding minimum energy
            fprintf(fp, "%f,%.8e,%.2f\n", a[i], Ealpha[i], w);
        }
        printf("Acceptance percentage is %f\n", 100.0*acpt/((double)N_al*(N_r-B)));
        printf("Minimum E(alpha) is: %f\nCorresponding to alpha: %f\n\n", Ealpha[idx], a[idx]);
    }
    fclose(fp);
    
    return 0;
}