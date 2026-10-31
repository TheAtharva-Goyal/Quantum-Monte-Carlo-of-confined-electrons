//=====================================================================================
// USING METROPOLIS-HASTINGS ALGORITHM TO CALCULATE MINIMUM ENERGY FOR TRIAL FUNCTION 1
//=====================================================================================


#include <stdio.h>
#include <math.h>
#include <stdint.h>
#include <time.h>

#define N_al 300 //number of steps for alpha
#define N_r 250000 //number of steps for the coordinates of electrons

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

double Psi(double r1, double r2, double w, double a){ //given wavefunction
    return exp(-a*w*(r1*r1 + r2*r2)/2);
}

double E(double r1, double r2, double a, double w){ //given local energy function
    return (w*w/2)*(r1*r1+r2*r2)*(1-(a*a))+3*a*w;
}

double radius(double x, double y, double z){
    return sqrt(x*x+y*y+z*z);
}

static double w=1.0; double a[N_al]={0.0}, Ealpha[N_al]={0.0};

int main(){
    rng_seed((uint64_t)time(NULL) ^ ((uint64_t)clock() << 32)); //initialising the RNG

    double step=1.44, stepa=0.4/(N_al-1); long int ar=0; //adjust step by hand to get around 50% acceptance
    for(int i=0; i<N_al; i++){
        a[i]=0.8+i*stepa;  //incrementing to the next alpha
        double Eloc=0.0; //initialising local energy for a given alpha, for a single coonfiguration of electrons
        //randomly picking initial values symmetrically near the origin:
        double cx1=rng_uniform()-0.5, cy1=rng_uniform()-0.5, cz1=rng_uniform()-0.5, cx2=rng_uniform()-0.5, cy2=rng_uniform()-0.5, cz2=rng_uniform()-0.5;
        cx1/=10; cx2/=10; cy1/=10; cy2/=10; cz1/=10; cz2/=10; //further narrowing the initial positions closer to the origin
        double cr1=radius(cx1, cy1, cz1), cr2=radius(cx2, cy2, cz2);
        for(int j=0; j<N_r; j++){
            double x1n=cx1+step*(rng_uniform()-0.5); //taking steps uniformly
            double y1n=cy1+step*(rng_uniform()-0.5);
            double z1n=cz1+step*(rng_uniform()-0.5);
            double x2n=cx2+step*(rng_uniform()-0.5);
            double y2n=cy2+step*(rng_uniform()-0.5);
            double z2n=cz2+step*(rng_uniform()-0.5);
            double r1n=radius(x1n, y1n, z1n), r2n=radius(x2n, y2n, z2n); 
            if(rng_uniform()<exp(-a[i]*w*(r1n*r1n + r2n*r2n - cr1*cr1 - cr2*cr2))){ //condition for moving forward
                cx1=x1n; cx2=x2n; cy1=y1n; cy2=y2n; cz1=z1n; cz2=z2n; //moving forward
                cr1=radius(cx1, cy1, cz1); cr2=radius(cx2, cy2, cz2);
                if(j>=N_r/5) ar++; //N_r/5 is the burn period
            }
            if(j<N_r/5) continue;
            Eloc+=E(cr1, cr2, a[i], w);
        }
        Ealpha[i]=1.25*Eloc/N_r; //storing the value of energy for the i'th alpha
        //E(alpha) is the expectation value of Eloc, where |Psi(r1, r2)|^2 is the probability
    }
    int idx=0;
    for(int i=0; i<N_al; i++){
        if(Ealpha[i]<Ealpha[idx]) idx=i; //finding the minimum of E(alpha)
    }
    printf("Acceptance percentage is %f\n", 125.0*ar/(double)(N_al*N_r));
    printf("Minimum E(alpha) is: %f\nCorresponding to alpha: %f\n", Ealpha[idx], a[idx]);
    return 0;
}