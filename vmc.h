#ifndef VMC_H
#define VMC_H
#include <stdint.h>
#include <math.h>

void rng_seed(uint64_t seed);
double rng_uniform(void);

typedef struct{double x1, y1, z1, x2, y2, z2, r1, r2, r12;} MH; //testing
void MH_init(MH *w); //testing
int MH_stepT1(MH *c, double step, double a, double w); //testing
int MH_stepT2(MH *c, double step, double a, double b, double w); //testing

double stepfinderT1(double a, double w, long long pilot); //testing
double stepfinderT2(double a, double b, double w, long long pilot); //testing

static inline double radius2(double x, double y, double z){return x*x+y*y+z*z;} //testing
static inline double logPT1(double r1, double r2, double w, double a){return -a*w*(r1+r2);} //testing
static inline double logPT2(double r1, double r2, double r12, double w, double a, double b){return -a*w*(r1+r2) + r12/(1+b*r12);} //testing
static inline double EL1(double r1, double r2, double r12, double a, double w){return (w*w/2)*(r1+r2)*(1-(a*a))+3*a*w+ 1/r12;} //testing
static inline double EL2(double r1, double r2, double r12, double a, double b, double w){ //testing
    double G=1+b*r12, K=0.5/(G*G);
    return EL1(r1, r2, r12, a, w)+K*((a*w*r12)-K-(2/r12)+(2*b/G));
}

#endif //VMC_H