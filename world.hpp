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
     * parameter file example
     * \code
     * delta_t 0.1
     * t_end 1.0
     * \endcode
     *
     * @param filename filename of the parameter file
     */
    void read_Parameter(const std::string &filename);

    /**
     * @brief read the particles from the given data file
     *
     * @param filename filename of the particle data file
     */
    void read_Particles(const std::string &filename);

    // unknown marks, that there is no treatment of the boarder, leaving indicates, that particles can escape of our world
    enum BorderType { unknown = 0, leaving = 1 }; 
    
    // Value-Defintions of the different String values
    enum Option { NAME, DELTA_T, T_END, LENGTH, UPPER_BORDER, LOWER_BORDER}; 
    
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
