//=====================================================================================
// USING METROPOLIS-HASTINGS ALGORITHM TO CALCULATE MINIMUM ENERGY FOR TRIAL FUNCTION 1
//=====================================================================================


#include <stdio.h>
#include <math.h>
#include <stdint.h>
#include <time.h>
#include "vmc.h"

#define N_al 201 //number of steps for alpha
#define N_r 300000 //number of steps for the coordinates of electrons
#define B (N_r/5) //burn-in

static inline double E(double r1, double r2, double a, double w){
    return (w*w/2)*(r1+r2)*(1-(a*a))+3*a*w;
}

static double w=1.0; double a[N_al]={0.0}, Ealpha[N_al]={0.0};

int main(){
    rng_seed((uint64_t)time(NULL) ^ ((uint64_t)clock() << 32)); //initialising the RNG

    double stepa=0.4/(N_al-1); long int ar=0;
    for(int i=0; i<N_al; i++){
        a[i]=0.8+i*stepa; //incrementing to the next alpha
        double Eloc=0.0; //initialising local energy for a given alpha, for a single coonfiguration of electrons
        double step=stepfinderT1(a[i], w, 10000);
        MH k; MH_init(&k);
        for(int j=0; j<N_r; j++){
            int ch=MH_stepT1(&k, step, a[i], w); if(j<B) continue; ar+=ch;
            Eloc+=E(k.r1, k.r2, a[i], w);
        }
        Ealpha[i]=1.25*Eloc/N_r; //storing the value of energy for the i'th alpha
        //E(alpha) is the expectation value of Eloc, where |Psi(r1, r2)|^2 is the probability
    }
    int idx=0;
    for(int i=0; i<N_al; i++){
        if(Ealpha[i]<Ealpha[idx]) idx=i; //finding the minimum of E(alpha)
    }
    printf("Acceptance percentage is %f\n", 100.0*ar/(double)(N_al*(N_r-B)));
    printf("Minimum E(alpha) is: %f\nCorresponding to alpha: %f\n", Ealpha[idx], a[idx]);
    return 0;
}