#include "vmc.h"

//---------xoshiro256** algorithm--------------------
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
//---------------------------------------------------

void MH_init(MH *w){
    w->x1=(rng_uniform()-0.5)/4; w->y1=(rng_uniform()-0.5)/4; w->z1=(rng_uniform()-0.5)/4;
    w->x2=(rng_uniform()-0.5)/4; w->y2=(rng_uniform()-0.5)/4; w->z2=(rng_uniform()-0.5)/4;
    w->r1=radius2(w->x1, w->y1, w->z1);
    w->r2=radius2(w->x2, w->y2, w->z2);
    w->r12=sqrt(radius2(w->x1-w->x2, w->y1-w->y2, w->z1-w->z2));
}

int MH_stepT1(MH *c, double step, double a, double w){
    double x1n=c->x1+step*(rng_uniform()-0.5), y1n=c->y1+step*(rng_uniform()-0.5), z1n=c->z1+step*(rng_uniform()-0.5);
    double x2n=c->x2+step*(rng_uniform()-0.5), y2n=c->y2+step*(rng_uniform()-0.5), z2n=c->z2+step*(rng_uniform()-0.5);
    double r1n=radius2(x1n, y1n, z1n), r2n=radius2(x2n, y2n, z2n), rn=sqrt(radius2(x1n-x2n, y1n-y2n, z1n-z2n));
    if(log(rng_uniform())<logPT1(r1n, r2n, w, a)-logPT1(c->r1, c->r2, w, a)){
        c->x1=x1n; c->y1=y1n; c->z1=z1n; c->x2=x2n; c->y2=y2n; c->z2=z2n;
        c->r1=r1n; c->r2=r2n; c->r12=rn;
        return 1;
    }
    return 0;
}

int MH_stepT2(MH *c, double step, double a, double b, double w){
    double x1n=c->x1+step*(rng_uniform()-0.5), y1n=c->y1+step*(rng_uniform()-0.5), z1n=c->z1+step*(rng_uniform()-0.5);
    double x2n=c->x2+step*(rng_uniform()-0.5), y2n=c->y2+step*(rng_uniform()-0.5), z2n=c->z2+step*(rng_uniform()-0.5);
    double r1n=radius2(x1n, y1n, z1n), r2n=radius2(x2n, y2n, z2n), rn=sqrt(radius2(x1n-x2n, y1n-y2n, z1n-z2n));
    if(log(rng_uniform())<logPT2(r1n, r2n, rn, w, a, b)-logPT2(c->r1, c->r2, c->r12, w, a, b)){
        c->x1=x1n; c->y1=y1n; c->z1=z1n; c->x2=x2n; c->y2=y2n; c->z2=z2n;
        c->r1=r1n; c->r2=r2n; c->r12=rn;
        return 1;
    }
    return 0;
}

double stepfinderT1(double a, double w, long long pilot){
    long long burn=pilot/2, ct=pilot-burn; //number of iterations of the loop and burn-in
    MH c; MH_init(&c); double delmin=0.0, delmax=10.0/sqrt(a*w), step;
    for(int i=0; i<10; i++){ //runs a maximum of 10 loops for tuning step size
        step=(delmax+delmin)/2.0; long long acpt=0;
        for(int j=0; j<pilot; j++){ //Metro-Hast loops
            int chk=MH_stepT1(&c, step, a, w); if(j>=burn) acpt+=chk;
        }
        double ar=100.0*acpt/ct;
        if(fabs(ar-50)<1.0) break;
        if(ar>50.0) delmin=step; else delmax=step;
    }
    return step;
}

double stepfinderT2(double a, double b, double w, long long pilot){
    long long burn=pilot/2, ct=pilot-burn; //number of iterations of the loop and burn-in
    MH c; MH_init(&c); double delmin=0.0, delmax=10.0/sqrt(a*w), step;
    for(int i=0; i<10; i++){ //runs a maximum of 10 loops for tuning step size
        step=(delmax+delmin)/2.0; long long acpt=0;
        for(int j=0; j<pilot; j++){ //Metro-Hast loops
            int chk=MH_stepT2(&c, step, a, b, w); if(j>=burn) acpt+=chk;
        }
        double ar=100.0*acpt/ct;
        if(fabs(ar-50)<1.0) break;
        if(ar>50.0) delmin=step; else delmax=step;
    }
    return step;
}