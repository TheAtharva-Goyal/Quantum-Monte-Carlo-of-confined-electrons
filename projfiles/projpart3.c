//========================================================================================================
// USING METROPOLIS-HASTINGS TO FIND MEAN SEPARATION OF ELECTRONS AND VARIANCE OF ENERGY AT ENERGY MINIMUM
//========================================================================================================


#include <stdio.h>
#include <math.h>
#include <stdint.h>
#include <time.h>
#include <stdlib.h>

#define N_al 500 //number of steps for alpha
#define N_r 2000000 //number of steps for the coordinates
#define B (N_r/5) //burn-in

//-----------xoshiro256** algo:-----------------
static uint64_t s[4];

static uint64_t splitmix64(uint64_t *x) {
    uint64_t z = (*x += 0x9E3779B97F4A7C15ULL);
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
    return z ^ (z >> 31);
}

static inline uint64_t rotl(uint64_t x, int k) {
    return (x << k) | (x >> (64 - k));
}

static uint64_t next_u64(void) {
    uint64_t result = rotl(s[1] * 5, 7) * 9;
    uint64_t t = s[1] << 17;
    s[2] ^= s[0]; s[3] ^= s[1];
    s[1] ^= s[2]; s[0] ^= s[3];
    s[2] ^= t;
    s[3] = rotl(s[3], 45);
    return result;
}

void rng_seed(uint64_t seed) {
    for (int i = 0; i < 4; i++) s[i] = splitmix64(&seed);
}

double rng_uniform(void) { //uniform in (0,1), 53-bit
    return ((next_u64() >> 11) + 0.5) * (1.0 / 9007199254740992.0);
}
//---------------------------------------

double PT(double r1, double r2, double w, double a){ //given wavefunction
    return exp(-a*w*(r1+r2)/2); //directly taking r1 and r2 as squares for computational ease (since radius2() outputs squares)
}

double EL(double r1, double r2, double a, double w, double r12){ //given local energy function
    return (w*w/2)*(r1+r2)*(1-(a*a))+3*a*w+ 1/r12;
}

double radius2(double x, double y, double z){ //square of (x,y,z) from origin
    return x*x+y*y+z*z;
}

static double w=1.0; double a[N_al]={0.0}, Ealpha[N_al]={0.0}, VarEalpha[N_al]={0.0}; //note: manually change w to 0.01 and 0.5 and run the programme

double stepfinder(double a, double w, double delmin, double delmax){ //finds the step size for ~50% acceptance for any pair (a,w)
    long long plt=B, burn=0.5*B, ct=plt-burn; //defining the number of iterations of the loop and burn-in
    double cx1=rng_uniform()-0.5, cy1=rng_uniform()-0.5, cz1=rng_uniform()-0.5; //randomly picking initial values symmetrically near the origin
    double cx2=rng_uniform()-0.5, cy2=rng_uniform()-0.5, cz2=rng_uniform()-0.5;
    double cr1=radius2(cx1, cy1, cz1), cr2=radius2(cx2, cy2, cz2);
    short ctr=1; double temp=100, step;
    while(fabs(temp-50.0)>0.1 && ctr<1000){
        long long acpt=0; step=(delmax+delmin)/2;
        for(int j=0; j<plt; j++){ //pilot run for tuning step size according to alpha
            double x1n=cx1+step*(rng_uniform()-0.5); //taking steps uniformly
            double y1n=cy1+step*(rng_uniform()-0.5);
            double z1n=cz1+step*(rng_uniform()-0.5);
            double x2n=cx2+step*(rng_uniform()-0.5);
            double y2n=cy2+step*(rng_uniform()-0.5);
            double z2n=cz2+step*(rng_uniform()-0.5);
            double r1n=radius2(x1n, y1n, z1n), r2n=radius2(x2n, y2n, z2n);
            if(log(rng_uniform())<(-a*w*(r1n+r2n-cr1-cr2))){
                cx1=x1n; cx2=x2n; cy1=y1n; cy2=y2n; cz1=z1n; cz2=z2n;
                cr1=radius2(cx1, cy1, cz1); cr2=radius2(cx2, cy2, cz2);
                if(j>=burn) acpt++;
            }
        }
        temp=100.0*acpt/ct;
        if(temp>50.0) delmin=step;
        if(temp<50.0) delmax=step;
        ctr++;
    }
    return step;
}

int main(){
    rng_seed((uint64_t)time(NULL) ^ ((uint64_t)clock() << 32)); //initialising the RNG
    
    double stepa=0.7/(N_al-1); long long acpt=0, ct=0;
    for(int i=0; i<N_al; i++){
        a[i]=0.3+i*stepa;  //incrementing to the next alpha with fixed interval
        double Eloc=0.0, E2=0.0; //total local energy for a given alpha
        double cx1=rng_uniform()-0.5, cy1=rng_uniform()-0.5, cz1=rng_uniform()-0.5; //randomly picking initial values symmetrically near the origin
        double cx2=rng_uniform()-0.5, cy2=rng_uniform()-0.5, cz2=rng_uniform()-0.5;
        double cr1=radius2(cx1, cy1, cz1), cr2=radius2(cx2, cy2, cz2);
        double step=stepfinder(a[i], w, 0.0, 10.0/sqrt(a[i]*w));
        for(int j=0; j<N_r; j++){
            double x1n=cx1+step*(rng_uniform()-0.5); //taking steps uniformly
            double y1n=cy1+step*(rng_uniform()-0.5);
            double z1n=cz1+step*(rng_uniform()-0.5);
            double x2n=cx2+step*(rng_uniform()-0.5);
            double y2n=cy2+step*(rng_uniform()-0.5);
            double z2n=cz2+step*(rng_uniform()-0.5);
            double r1n=radius2(x1n, y1n, z1n), r2n=radius2(x2n, y2n, z2n); 
            if(log(rng_uniform())<(-a[i]*w*(r1n+r2n-cr1-cr2))){ //condition for moving forward
                cx1=x1n; cx2=x2n; cy1=y1n; cy2=y2n; cz1=z1n; cz2=z2n;
                cr1=r1n; cr2=r2n;
                if(j>=B) acpt++;
            }
            if(j<B) continue; double r12=radius2(cx1-cx2, cy1-cy2, cz1-cz2); if(r12==0.0) continue; //(since EL has a 1/r12 term)
            double e=EL(cr1, cr2, a[i], w, sqrt(r12)); Eloc+=e; E2+=e*e;
        }
        Ealpha[i]=Eloc/(N_r-B); //storing the value of E(alpha)
        VarEalpha[i]=E2/(N_r-B)-(Ealpha[i]*Ealpha[i]); //variance of E(alpha)
    }

    int idx=0;
    for(int i=0; i<N_al; i++){
        if(Ealpha[i]<Ealpha[idx]) idx=i;
    }
    printf("Acceptance percentage is %f\n", 100.0*acpt/(N_al*(N_r-B)));
    printf("Minimum E(alpha) is: %f\nCorresponding to alpha: %f\n", Ealpha[idx], a[idx]);

    //Computing expectation value of mean distance between electrons at energy minimum:
    double al=a[idx], r=0.0, step=stepfinder(al, w, 0.0, 10.0/sqrt(al*w));
    double cx1=rng_uniform()-0.5, cy1=rng_uniform()-0.5, cz1=rng_uniform()-0.5;
    double cx2=rng_uniform()-0.5, cy2=rng_uniform()-0.5, cz2=rng_uniform()-0.5;
    double cr1=radius2(cx1, cy1, cz1), cr2=radius2(cx2, cy2, cz2);
    for(int j=0; j<N_r; j++){
        double x1n=cx1+step*(rng_uniform()-0.5);
        double y1n=cy1+step*(rng_uniform()-0.5);
        double z1n=cz1+step*(rng_uniform()-0.5);
        double x2n=cx2+step*(rng_uniform()-0.5);
        double y2n=cy2+step*(rng_uniform()-0.5);
        double z2n=cz2+step*(rng_uniform()-0.5);
        double r1n=radius2(x1n, y1n, z1n), r2n=radius2(x2n, y2n, z2n); 
        if(log(rng_uniform())<(-al*w*(r1n+r2n-cr1-cr2))){
            cx1=x1n; cx2=x2n; cy1=y1n; cy2=y2n; cz1=z1n; cz2=z2n;
            cr1=r1n; cr2=r2n;
        }
        if(j<B) continue;
        r+=sqrt(radius2(cx1-cx2, cy1-cy2, cz1-cz2));
    }
    r/=(N_r-B);
    printf("Expectation of mean separation for\nw=%.2f and alpha=%f is:\n<r12> = %f\nVar(E) at the minimum is: %f\n", w, al, r, VarEalpha[idx]);

    //plotting the graps using Python: 
    FILE *fp=fopen("data.csv", "w"); if(!fp) return 1;
    for(int i=0; i<N_al; i++) fprintf(fp, "%f,%f,%f\n", a[i], Ealpha[i], VarEalpha[i]);
    fclose(fp);
    char cmd[100];
    sprintf(cmd, "uv run plotter.py %.2f", w); //command to run the python file taking input value of w
    system(cmd);

    return 0;
}