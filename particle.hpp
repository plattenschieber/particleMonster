#ifndef _PARTICLE_HPP
#define _PARTICLE_HPP

#include "defines.hpp"

/**
 * @brief Particle data
 *
 * this class contains the particle data
 */
class Particle {
public:
    /// ID
    int ID;
    /// Mass
    real m;
    /// Position
    real x[DIM];
    /// Velocity
    real v[DIM];
    /// Force
    real F[DIM];
    /// Force (previous step)
    real F_old[DIM];
	/// Clear all data of the particle (for tmp particles)
	void clear() {
		m = 0;
		for(int i=0; i<DIM; i++)
		{
			x[i] = 0;
			v[i] = 0;
			F[i] = 0;
			F_old[i] = 0;
		}
	}
};

#endif // _PARTICLE_HPP
// vim:set et sts=4 ts=4 sw=4 ai ci cin:
