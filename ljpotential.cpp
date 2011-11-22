#include "ljpotential.hpp"
#include <math.h>
#include <iostream>

real LJPotential::force(Particle &p, Particle &q) /*, real eps, real sigma*/
{
    real eps = 1.0, sigma = 1.0;
    real rijsq = 0.0;
    for (unsigned int d=0; d<DIM; d++)
        rijsq += sqr(q.x[d]-p.x[d]);
    
    
    // DEBUG std::cout << "rijsq: " << rijsq << "###### pow: " << pow(sigma*sigma/rijsq,3) << std::endl;
    real Urij = 4*eps*pow(sigma*sigma/rijsq,3)*(pow(sigma*sigma/rijsq,3) - 1);
    
	for (int d=0; d<DIM; d++) {
		double force = 24*eps*(1/rijsq)*pow(sigma*sigma/rijsq,3)*(1-2*pow(sigma*sigma/rijsq,3))*(q.x[d]-p.x[d]);
        p.F[d] += force;
        q.F[d] -= force;
    }
	return Urij;
}

// vim:set et sts=4 ts=4 sw=4 ai ci cin:
