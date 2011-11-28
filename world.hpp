#ifndef _WORLD_HPP
#define _WORLD_HPP

#include "defines.hpp"
#include "particle.hpp"
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <map>

/**
 * @brief the world class holds all information of the simulation environment
 */


class World {
public:
    World();

    /**
     * @brief read the world parameters from the given parameter file
     *
     * parameter file example TODO: rewrite code example
     * \code
     * delta_t 0.1
     * t_end 1.0
     * \endcode
     *
     * @param filename filename of the parameter file
     */
    virtual void read_Parameter(const std::string &filename);

    /**
     * @brief read the particles from the given data file
     *
     * @param filename filename of the particle data file
     */
    virtual void read_Particles(const std::string &filename);

    // unknown marks, that there is no treatment of the boarder, leaving indicates, that particles can escape of our world and periodic will let the particles enter on the opposite side
    enum BorderType { unknown = 0, leaving = 1, periodic = 2 }; 
    
    // Value-Defintions of the different String values
    enum Option { NAME, DELTA_T, T_END, LENGTH, UPPER_BORDER, LOWER_BORDER, EPSILON, SIGMA}; 
    
    // Map to associate the strings with the enum values
    std::map<std::string, World::Option> mapOptions;
    // data structures
    /// Name of the simulated world
    std::string name;
    /// Current time
    real t;
    /// Timestep
    real delta_t;
    /// End of simulation
    real t_end;
    /// kinetic energy
    real e_kin;
    /// potential energy
    real e_pot;
    /// total energy
    real e_tot;
    /// the axis lengths of our world
    real length[DIM];
    /// zero breakthrough
    real sigma;
    // potential depth
    real epsilon;
    /// Vector of particles
    std::vector<Particle> particles;
    /// upper borders 
    BorderType upper_border[DIM];
    /// lower borders
    BorderType lower_border[DIM];
};

/**
 * @brief a ostream operator for the World class
 *
 * @param os stream object
 * @param W the world
 *
 * @return resulting stream object
 */
std::ostream& operator << (std::ostream& os, World& W);

#endif // _WORLD_HPP
// vim:set et sts=4 ts=4 sw=4 ai ci cin:
