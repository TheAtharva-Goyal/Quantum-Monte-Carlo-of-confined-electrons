//=======================================================================================================
// USING METRO-HAST ALGORITHM TO PLOT ENERGY VS ALPHA AT OPTIMAL BETA AND ENERGY VS BETA AT OPTIMAL ALPHA
// FOR TRIAL FUNCTION 2
//=======================================================================================================


#include <stdio.h>
#include <math.h>
#include <stdint.h>
#include <time.h>
#include "vmc.h"

#define N_r 200000 //number of steps for the coordinates
#define N_f 50000000 //number of steps for final evaluation
#define B (N_r/4) //burn-in
#define PLT 40000 //number of loops for stepfinder
#define N_l 11 //number of values checked for alpha and beta
#define N_plot 301 //number of values of alpha and beta for graph plotting
#define Aspan 0.10 //width of range of alpha and beta for graph plotting
#define Bspan 0.10

static double wval[3]={0.01, 0.5, 1.0}; double a0val[3]={0.438878, 0.824337, 0.882565}, vals[3]={0.0}; //vals stores the results of E, Var(E), <r12>
//wval are the required values of w
//a0val are the optimal values of alpha from the previous part for the respective w's

void eval_step(double a, double b, double w, long long N, double step){ //calculates E, Var(E), <r12>
    double En=0.0, En2=0.0, mr=0.0; int bn=N/4;
    MH k; MH_init(&k);
    for(int j=0; j<N; j++){
        MH_stepT2(&k, step, a, b, w); if(j<bn) continue;
        double e=EL2(k.r1, k.r2, k.r12, a, b, w);
        En+=e; En2+=e*e; mr+=k.r12;
    }
    vals[0]=En/(N-bn); vals[1]=En2/(N-bn)-vals[0]*vals[0]; vals[2]=mr/(N-bn);
}

void eval(double a, double b, double w, long long N){ //calculates E, Var(E), <r12>
    double step=stepfinderT2(a, b, w, PLT);
    eval_step(a, b, w, N, step);
}

int main(){
    rng_seed((uint64_t)time(NULL) ^ ((uint64_t)clock() << 32)); //initialising the RNG

    FILE *alpha_file=fopen("energy_vs_alpha.csv", "w");
    FILE *beta_file=fopen("energy_vs_beta.csv", "w");
    if(alpha_file==NULL || beta_file==NULL){
        perror("Error opening CSV file");
        if(alpha_file!=NULL) fclose(alpha_file);
        if(beta_file!=NULL) fclose(beta_file);
        return 1;
    }
    
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
                double ca=a-da+k*sa; if(ca<0.0) continue;
                eval(ca, b, w, N_r);
                if(vals[0]<best){best=vals[0]; ab=ca; ia=k;}
            }
            a=ab; //taking the best found value of alpha

            if(ia>0 && ia<N_l-1) da/=2.0; //narrowing the range of alpha and beta if minimum found
            if(ib>0 && ib<N_l-1) db/=2.0;
        }
        eval(a, b, w, N_f);
        printf("for w=%.2f:\nalpha=%f\nbeta=%f\nE=%f\nvarE=%f\nmeanr=%f\n\n", w, a, b, vals[0], vals[1], vals[2]);

        double pltst=stepfinderT2(a, b, w, 2*PLT);
        int mid=N_plot/2; double alpha_down=fmin(Aspan,a);
        for(int x=0; x<N_plot; x++){
            double ca;
            if(x<=mid) ca=a-alpha_down*(double)(mid-x)/mid;
            else ca=a+Aspan*(double)(x-mid)/mid;
            eval_step(ca, b, w, N_r, pltst); fprintf(alpha_file, "%.2f,%.8e,%.8e,%.8e,%.12e\n", w, ca, b, a, vals[0]);
        }

        double beta_down=fmin(Bspan, b);
        for(int x=0; x<N_plot; x++){
            double cb;
            if(x<=mid) cb=b-beta_down*(double)(mid-x)/mid;
            else cb=b+Bspan*(double)(x-mid)/mid;
            eval_step(a, cb, w, 5*N_r, pltst); fprintf(beta_file, "%.2f,%.8e,%.8e,%.8e,%.12e\n", w, cb, a, b, vals[0]);
        }
    }
    fclose(alpha_file); fclose(beta_file);

    return 0;
}