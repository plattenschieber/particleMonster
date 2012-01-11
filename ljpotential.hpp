#ifndef _LJPOTENTIAL_HPP
#define _LJPOTENTIAL_HPP

#include "potential.hpp"

/**
 * @brief A pseudo Potential called Lennard-Jones        
 */
class LJPotential : public Potential {
public:
    /**
     * @brief calculate the force bewteen the two particles and add it to p
     *
     * @param p particle p
     * @param q particle q
     * @param eps depth of potential
     * @param sigma equals the zero paththrough
     *
     * @return potential energy
     */
    real force (Particle& p, Particle& q, real dist, real distV[], real eps, real sigma);

private:
    real force (Particle& p, Particle& q){return 0.0;};
};

#endif // _LJPOTENTIAL_HPP
// vim:set et sts=4 ts=4 sw=4 ai ci cin:
