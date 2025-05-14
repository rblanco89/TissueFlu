/*=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- RANDOM NUMBERS STRUCTURES =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
// These generators were taken from "Numerical Recipes: The Art of Scientific Computing, 3rd Ed"
#include <vector>

/*===============================================================================================*/
// UNIFORM
/*===============================================================================================*/
// This is a recommend uniform random generator, its period is ~ 3.138 x 10^57.
struct Ran
{
        unsigned long long u,v,w;
        // Call with any integer seed (exept value of v below)
        Ran(unsigned long long j) : v(4101842887655102017LL), w(1)
        {
                u = j^v; int64();
                v = u; int64();
                w = v; int64();
        }
        // Return 64-bit random integer
        inline unsigned long long int64()
        {
                u = u * 2862933555777941757LL + 7046029254386353087LL;
                v ^= v >> 17; v ^= v << 31; v ^= v >> 8;
                w = 4294957665U*(w & 0xffffffff) + (w >> 32);
                unsigned long long x = u ^ (u << 21); x ^= x >> 35; x ^= x << 4;
                return (x + v) ^ w;
        }
        // Return random double-precision floating value in the range from 0 to 1
        inline double doub() { return 5.42101086242752217E-20 * int64(); } // multiply by 1/INT_MAX
        // Return 32-bit random integer
        inline unsigned int int32() { return (unsigned int)int64(); }
};

/*===============================================================================================*/
// NORMAL 
/*===============================================================================================*/

struct Normaldev : Ran {
        double mu,sig;
        Normaldev(double mmu, double ssig, unsigned long long i) : Ran(i), mu(mmu), sig(ssig){}
        double dev() {
                double u,v,x,y,q;
                do {
                        u = doub();
                        v = 1.7156*(doub()-0.5);
                        x = u - 0.449871;
                        y = abs(v) + 0.386595;
                        q = x*x + y*(0.19600*y-0.25472*x);
                } while (q > 0.27597 && (q > 0.27846 || v*v > -4.*log(u)*u*u));
                return mu + sig*v/u;
        }
};

/*===============================================================================================*/
// POISSON
/*===============================================================================================*/

double gammln(double xx) {
	int j;
	double x, y, tmp, ser;
	static const double cof[14]={57.1562356658629235,-59.5979603554754912,
	14.1360979747417471,-0.491913816097620199,.339946499848118887e-4,
	.465236289270485756e-4,-.983744753048795646e-4,.158088703224912494e-3,
	-.210264441724104883e-3,.217439618115212643e-3,-.164318106536763890e-3,
	.844182239838527433e-4,-.261908384015814087e-4,.368991826595316234e-5};
	
	if (xx <= 0) throw("bad arg in gammln");

	y = x = xx;
	tmp = x + 5.24218750000000000; // Rational 671/128.
	tmp = (x + 0.5)*log(tmp) - tmp;
	ser = 0.999999999999997092;
	for (j=0; j<14; j++) ser += cof[j]/++y;
	return tmp + log(2.5066282746310005*ser/x);
}

struct Poissondev : Ran {
	double lambda, sqlam, loglam, lamexp, lambold;
	std::vector<double> logfact;
	int swch;
	Poissondev(double llambda, unsigned long long i) : Ran(i), lambda(llambda),
	logfact(1024,-1.), lambold(-1.) {}
	int dev() {
		// Return a Poisson deviate using the most recently set value of lambda
		double u, u2, v, v2, p, t, lfac;
		int k;
		if (lambda < 5.) {
			//Will use product of uniforms method.
			if (lambda != lambold) lamexp = exp(-lambda);
			k = -1;
			t = 1.0;
			do {
				++k;
				t *= doub();
			} while (t > lamexp);
		} else {
			// Will use ratio-of-uniforms method.
			if (lambda != lambold) {
				sqlam = sqrt(lambda);
				loglam = log(lambda);
			}
			for (;;) {
				u = 0.64*doub();
				v = -0.68 + 1.28*doub();
				if (lambda > 13.5) { // Outer squeeze for fast rejection.
					v2 = v*v;
					if (v >= 0.0) {if (v2 > 6.5*u*(0.64-u)*(u+0.2)) continue;}
					else {if (v2 > 9.6*u*(0.66-u)*(u+0.07)) continue;}
				}
				k = int(floor(sqlam*(v/u)+lambda+0.5));
				if (k < 0) continue;
				u2 = u*u;
				if (lambda > 13.5) { // Inner squeeze for fast acceptance.
					if (v >= 0.) {if (v2 < 15.2*u2*(0.61-u)*(0.8-u)) break;}
					else {if (v2 < 6.76*u2*(0.62-u)*(1.4-u)) break;}
				}
				if (k < 1024) {
					if (logfact[k] < 0.) logfact[k] = gammln(k+1.);
					lfac = logfact[k];
				} else lfac = gammln(k+1.);
				p = sqlam*exp(-lambda + k*loglam - lfac); // Only when we must.
				if (u2 < p) break;
			}
		}
		lambold = lambda;
		return k;
	}
	int dev(double llambda) {
		// Reset lambda and then return a Poisson deviate.
		lambda = llambda;
		return dev();
	}
};

/*===============================================================================================*/
// GAMMA 
/*===============================================================================================*/

struct Gammadev : Normaldev {
	double alph, oalph, bet;
	double a1, a2;
	Gammadev (double aalph, double bbet, unsigned long long i)
		: Normaldev(0.0,1.0,i), alph(aalph), oalph(aalph), bet(bbet) {
		if (alph <= 0.0) throw("Bad alpha in Gammadev");
		if (alph < 1.0) alph += 1.0;
		a1 = alph - 1.0/3.0;
		a2 = 1.0/sqrt(9.0*a1);
	}
	double dev() { 
		double u, v, x;
		do {
			do {
				x = Normaldev::dev();
				v = 1.0 + a2*x;
			} while (v <= 0.0);
			v = v*v*v;
			u = doub();
		} while (u > 1.0 - 0.331*x*x*x*x && log(u) > 0.5*x*x + a1*(1.0 - v + log(v)));
		if (alph == oalph) return a1*v/bet;
		else {
			do u = doub(); while (u == 0.0);
			return pow(u,1.0/oalph)*a1*v/bet;
		}
	}
};
