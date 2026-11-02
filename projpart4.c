//======================================================================================================
// USING METROPOLIS-HASTINGS ALGORITHM ON TRIAL FUNCTION 2 TO FIND MEAN ELECTRON SEPARATION AND VARIANCE
// AT ENERGY MINIMUM FOR OMEGA = 0.01, 0.5 and 1.0
//======================================================================================================


#include <stdio.h>
#include <math.h>
#include <stdint.h>
#include <time.h>
#include "vmc.h"

#define N_r 50000 //number of steps for the coordinates
#define N_f 5000000 //number of steps for final evaluation
#define B (N_r/4) //burn-in
#define PLT 30000 //number of loops for stepfinder
#define N_l 11 //number of values checked for alpha and beta

static double wval[3]={0.01, 0.5, 1.0}; double a0val[3]={0.438878, 0.824337, 0.882565}, vals[3]={0.0}; //vals stores the results of E, Var(E), <r12>
//wval are the required values of w
//a0val are the optimal values of alpha from the previous part for the respective w's

void eval(double a, double b, double w, long long N){ //calculates E, Var(E), <r12>
    double step=stepfinderT2(a, b, w, PLT), En=0.0, En2=0.0, mr=0.0; int bn=N/4;
    MH k; MH_init(&k);
    for(int j=0; j<N; j++){
        MH_stepT2(&k, step, a, b, w); if(j<bn) continue;
        double e=EL2(k.r1, k.r2, k.r12, a, b, w);
        En+=e; En2+=e*e; mr+=k.r12;
    }
    vals[0]=En/(N-bn); vals[1]=En2/(N-bn)-vals[0]*vals[0]; vals[2]=mr/(N-bn);
}

int main(){
    rng_seed((uint64_t)time(NULL) ^ ((uint64_t)clock() << 32)); //initialising the RNG
    
    for(int i=0; i<3; i++){ //one loop for each value of w
        double w=wval[i], a=a0val[i], b=1.0; //starting (middle) values of alpha and beta
        double da=0.4, db=1.0; //half of range of alpha and beta
        for(int j=0; j<8; j++){
            double sa=2*da/(N_l-1), sb=2*db/(N_l-1); //step size for alpha and beta
            double best=1.0e20, bb=b; int ib=0; //finding the best value of beta and its index
            for(int k=0; k<N_l; k++){
                double cb=b-db+k*sb; if(cb<0.0) continue;
                eval(a, cb, w, N_r);
                if(vals[0]<best){best=vals[0]; bb=cb; ib=k;}
            } 
            b=bb; //taking the best found value of beta

            best=1.0e20; double ab=a; int ia=0; //finding the best value of alpha and its index
            for(int k=0; k<N_l; k++){
                double ca=a-da+k*sa; if(ca<=0.0) continue;
                eval(ca, b, w, N_r);
                if(vals[0]<best){best=vals[0]; ab=ca; ia=k;}
            }
            a=ab; //taking the best found value of alpha

            if(ia>0 && ia<N_l-1) da/=2.0; //narrowing the range of alpha and beta if minimum found
            if(ib>0 && ib<N_l-1) db/=2.0;
        }
        eval(a, b, w, N_f);
        printf("for w=%.2f:\nalpha=%f\nbeta=%f\nE=%f\nvarE=%f\nmeanr=%f\n\n", w, a, b, vals[0], vals[1], vals[2]);
    }

    return 0;
}