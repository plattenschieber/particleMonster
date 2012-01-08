#ifndef _WORLD_HPP
#define _WORLD_HPP

#include "defines.hpp"
#include "particle.hpp"
#include <vector>
#include <list>
#include <string>
#include <iostream>
#include <fstream>
#include <map>
#include <stdexcept>
#include <sstream>
#include <cmath>

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
     * name example01
     * delta_t 0.1
     * t_end 1.0
     * epsilon 1
     * sigma 1
     * length 10 10 15
     * upper_border leaving leaving leaving
     * lower_border leaving leaving leaving
     * \endcode
     *
     * @param filename filename of the parameter file
     */
    virtual void readParameter(const std::string &filename);

    /**
     * @brief read the particles from the given data file
     *
     * @param filename filename of the particle data file
     */
    virtual void readParticles(const std::string &filename);

    /**
     * @brief calculate the new beta
     */
    virtual real calcBeta();


    // unknown marks, that there is no treatment of the boarder, leaving indicates, that particles can escape of our world and periodic will let the particles enter on the opposite side
    /// Type of World Border
    enum borderType { unknown = 0, leaving = 1, periodic = 2 };
    /// Value-Defintions of the different option strings
    // DEFAULT is needed to handle unknown options - otherwise a new option with value 0 is created and will map NAME 
    enum Option { DEFAULT=0, NAME=1, DELTA_T=2, T_END=3, LENGTH=4, UPPERBORDER=5, LOWERBORDER=6, EPSILON=7, SIGMA=8,
                  STARTTEMP=9, STEPINTERVAL=10, TARGETTEMP=11, RANDOMSEED=12 };
    /// Map to associate the strings with the enum values
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
    /// average total energy for max the last 100 values (if more are inserted, the first bumps out)
    std::list<real> e_avglist;
    /// the sum over the last 100 total energys
    real e_avg;
    /// the axis lengths of our world
    real length[DIM];
    /// zero breakthrough
    real sigma;
    /// potential depth
    real epsilon;
    /// Number of particles overall
    int nParticles;
    /// List of particles
    std::list<Particle> particles;
    /// upper borders 
    borderType upper_border[DIM];
    /// lower borders
    borderType lower_border[DIM];

    // Thermostat
    /// if set, the start velocity of all particles is set to it
    real thermo_start_temp;
    ///
    bool isThermoStartTemp;
    /// every
    real thermo_step_interval;
    real thermo_target_temp;

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
