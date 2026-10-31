//=======================================================================================================
// USING METRO-HAST ALGORITHM TO PLOT ENERGY VS ALPHA AT OPTIMAL BETA AND ENERGY VS BETA AT OPTIMAL ALPHA
// FOR TRIAL FUNCTION 2
//=======================================================================================================


#include <stdio.h>
#include <math.h>
#include <stdint.h>
#include <time.h>

#define N_r 500000 //number of steps for the coordinates
#define N_f 50000000 //number of steps for final evaluation
#define B (N_r/4) //burn-in
#define PLT 50000 //number of loops for stepfinder
#define N_l 21 //number of values checked for alpha and beta
#define N_plot 501 //number of values of alpha and beta for graph plotting
#define Aspan 0.10 //width of range of alpha and beta for graph plotting
#define Bspan 0.10

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

static inline double logPT(double r1, double r2, double r12, double w, double a, double b){ //log of abs val squared of given wavefunction
    return (-a*w*(r1+r2))+(r12/(1+b*r12)); //directly taking r1 and r2 as squares (since radius2 outputs squares)
}

static inline double EL(double r1, double r2, double r12, double a, double b, double w){ //given local energy function
    double EL1=(w*w/2)*(r1+r2)*(1-(a*a))+3*a*w+ 1/r12, G=1+b*r12, K=0.5/(G*G);
    return EL1+ K*((a*w*r12)-K-(2/r12)+(2*b/G)); //r1, r2 are squares, but r12 is not
}

static inline double radius2(double x, double y, double z){ //square of distance of (x,y,z) from origin
    return x*x+y*y+z*z;
}

typedef struct{double x1, y1, z1, x2, y2, z2, r1, r2, r12;}MH; //current values of the Metro-Hast walker

void MH_init(MH *w){ //initialises the Metro-Hast walker to values symmetrically near the origin
    w->x1=(rng_uniform()-0.5)/4; w->y1=(rng_uniform()-0.5)/4; w->z1=(rng_uniform()-0.5)/4;
    w->x2=(rng_uniform()-0.5)/4; w->y2=(rng_uniform()-0.5)/4; w->z2=(rng_uniform()-0.5)/4;
    w->r1=radius2(w->x1, w->y1, w->z1);
    w->r2=radius2(w->x2, w->y2, w->z2);
    w->r12=sqrt(radius2(w->x1-w->x2, w->y1-w->y2, w->z1-w->z2));
}

int MH_step(MH *c, double step, double a, double b, double w){ //step incrementing function for the Metro-Hast walker
    double x1n=c->x1+step*(rng_uniform()-0.5), y1n=c->y1+step*(rng_uniform()-0.5), z1n=c->z1+step*(rng_uniform()-0.5);
    double x2n=c->x2+step*(rng_uniform()-0.5), y2n=c->y2+step*(rng_uniform()-0.5), z2n=c->z2+step*(rng_uniform()-0.5);
    double r1n=radius2(x1n, y1n, z1n), r2n=radius2(x2n, y2n, z2n), rn=sqrt(radius2(x1n-x2n, y1n-y2n, z1n-z2n));
    if(log(rng_uniform())<logPT(r1n, r2n, rn, w, a, b)-logPT(c->r1, c->r2, c->r12, w, a, b)){
        c->x1=x1n; c->y1=y1n; c->z1=z1n; c->x2=x2n; c->y2=y2n; c->z2=z2n;
        c->r1=r1n; c->r2=r2n; c->r12=rn;
        return 1;
    }
    return 0;
}

static double wval[3]={0.01, 0.5, 1.0}; double a0val[3]={0.438878, 0.824337, 0.882565}, vals[3]={0.0}; //vals stores the results of E, Var(E), <r12>
//wval are the required values of w
//a0val are the optimal values of alpha from the previous part for the respective w's

double stepfinder(double a, double b, double w){ //finds the appropriate step size for ~50% acceptance for (a,b,w)
    long long burn=PLT/2, ct=PLT-burn; //number of iterations of the loop and burn-in
    MH c; MH_init(&c); double delmin=0.0, delmax=30.0, step;
    for(int i=0; i<10; i++){ //runs a maximum of 10 loops for tuning step size
        step=(delmax+delmin)/2.0; long long acpt=0;
        for(int j=0; j<PLT; j++){ //Metro-Hast loops
            int chk=MH_step(&c, step, a, b, w); if(j>=burn) acpt+=chk;
        }
        double ar=100.0*acpt/ct;
        if(fabs(ar-50)<0.5) break;
        if(ar>50.0) delmin=step; else delmax=step;
    }
    return step;
}

void eval_step(double a, double b, double w, long long N, double step){
    double En=0.0, En2=0.0, mr=0.0; int bn=N/4;
    MH k; MH_init(&k);
    for(int j=0; j<N; j++){
        MH_step(&k, step, a, b, w); if(j<bn) continue;
        double e=EL(k.r1, k.r2, k.r12, a, b, w);
        En+=e; En2+=e*e; mr+=k.r12;
    }
    vals[0]=En/(N-bn); vals[1]=En2/(N-bn)-vals[0]*vals[0]; vals[2]=mr/(N-bn);
}

void eval(double a, double b, double w, long long N){ //calculates E, Var(E), <r12>
    double step=stepfinder(a, b, w);
    eval_step(a, b, w, N, step);
}

int main(){
    rng_seed((uint64_t)time(NULL) ^ ((uint64_t)clock() << 32)); //initialising the RNG

    FILE *alpha_file=fopen("energy_vs_alpha.csv", "w");
    FILE *beta_file=fopen("energy_vs_beta.csv", "w");
    if (alpha_file==NULL || beta_file==NULL) {
        perror("Error opening CSV file");
        if(alpha_file!=NULL) fclose(alpha_file);
        if(beta_file!=NULL) fclose(beta_file);
        return 1;
    }
    fprintf(alpha_file, "omega,alpha,beta_fixed,alpha_final,E\n");
    fprintf(beta_file,  "omega,beta,alpha_fixed,beta_final,E\n");
    
    for(int i=0; i<3; i++){ //one loop for each value of w
        double w=wval[i], a=a0val[i], b=1.0; //starting (middle) values of alpha and beta
        double da=0.4, db=1.0; //half of range of alpha and beta
        for(int j=0; j<15; j++){
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

        double pltst=stepfinder(a, b, w);
        int mid=N_plot/2; double alpha_down=fmin(Aspan,a);
        for(int x=0; x<N_plot; x++){
            double ca;
            if(x<=mid) ca=a-alpha_down*(double)(mid-x)/mid;
            else ca=a+Aspan*(double)(x-mid)/mid;
            eval_step(ca, b, w, 20*N_r, pltst); fprintf(alpha_file, "%.8g,%.12g,%.12g,%.12g,%.12g\n", w, ca, b, a, vals[0]);
        }

        double beta_down=fmin(Bspan, b);
        for(int x=0; x<N_plot; x++){
            double cb;
            if(x<=mid) cb=b-beta_down*(double)(mid-x)/mid;
            else cb=b+Bspan*(double)(x-mid)/mid;
            eval_step(a, cb, w, 100*N_r, pltst); fprintf(beta_file, "%.8g,%.12g,%.12g,%.12g,%.12g\n", w, cb, a, b, vals[0]);
        }
    }
    fcloseall;

    return 0;
}