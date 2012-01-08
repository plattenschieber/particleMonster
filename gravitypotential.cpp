#include "gravitypotential.hpp"
#include <cmath>

real GravityPotential::force(Particle &p, Particle &q)
{
    // squared euclidian distance - r - between p and q
    real r2 = 0.0;
    for (unsigned int d=0; d<DIM; d++)
        r2 += sqr(q.x[d]-p.x[d]);
    // we don't need no sqrt to exponentiate r by 3 
    real r3 = r2 * sqrt(r2);
    // scalar factor of our directional force-vector (induced by scaled gravitypotential)
    real scal  = -p.m*q.m/r3;
    // add force to our particle p
    for (int d=0; d<DIM; d++)
	p.F[d] += -scal * (q.x[d]-p.x[d]);
    return scal;
}

// vim:set et sts=4 ts=4 sw=4 ai ci cin:
