#include "ljpotential.hpp"
#include <math.h>
#include <iostream>
real LJPotential::force(Particle &p, Particle &q, real dist , real eps, real sigma)
{
    // potential depth and zero breakthrough
    //real eps = 1.0, sigma = 1.0;
    
    // potential between p and q
    real Urij = 4*eps*pow(sqr(sigma)/dist,3)*(pow(sqr(sigma)/dist,3) - 1);
    
    // force is transfered from q to p by means of ljpotential
    for (int d=0; d<DIM; d++) {
	double force = 24*eps*(1/dist)*pow(sigma*sigma/dist,3)*(1-2*pow(sigma*sigma/dist,3))*(q.x[d]-p.x[d]);
	p.F[d] += force;
	q.F[d] -= force;
    }
    return Urij;
}

// vim:set et sts=4 ts=4 sw=4 ai ci cin:
