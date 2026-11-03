#ifndef VMC_H
#define VMC_H
#include <stdint.h>
#include <math.h>

void rng_seed(uint64_t seed); //part of the xoshiro256** RNG algorithm
double rng_uniform(void);

typedef struct{double x1, y1, z1, x2, y2, z2, r1, r2, r12;} MH; //current coordinates of the Metro-Hast walker
void MH_init(MH *w); //initialises the Metro-Hast walker to values symmetrically near the origin
int MH_stepT1(MH *c, double step, double a, double w); //step incrementing function for the Metro-Hast walker (corresponding to Psi_T1)
int MH_stepT2(MH *c, double step, double a, double b, double w); //step incrementing function for the Metro-Hast walker (corresponding to Psi_T2)

double stepfinderT1(double a, double w, long long pilot); //finds the appropriate step size for ~50% acceptance for (a,b)
double stepfinderT2(double a, double b, double w, long long pilot); //finds the appropriate step size for ~50% acceptance for (a,b,w)

static inline double radius2(double x, double y, double z){return x*x+y*y+z*z;} //square of distance of (x,y,z) from origin
static inline double logPT1(double r1, double r2, double w, double a){return -a*w*(r1+r2);} //log(|Psi_T1|^2)
static inline double logPT2(double r1, double r2, double r12, double w, double a, double b){return -a*w*(r1+r2) + r12/(1+b*r12);} //log(|Psi_T2|^2)
static inline double EL1(double r1, double r2, double r12, double a, double w){return (w*w/2)*(r1+r2)*(1-(a*a))+3*a*w+ 1/r12;} //local energy function corresponding to Psi_T1
static inline double EL2(double r1, double r2, double r12, double a, double b, double w){ //local energy function corresponding to Psi_T2
    double G=1+b*r12, K=0.5/(G*G);
    return EL1(r1, r2, r12, a, w)+K*((a*w*r12)-K-(2/r12)+(2*b/G));
}

#endif //VMC_H