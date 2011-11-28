#ifndef _WORLDLC_HPP
#define _WORLDLC_HPP

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
class WorldLC : public World {
public:
    WorldLC();

    /* @brief read the worldlc parameters from the given parameter file
     *
     * parameter file example TODO: rewrite code example
     * \code
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

    // Value-Defintions of the different String values
    enum Option { NAME, DELTA_T, T_END, LENGTH, UPPER_BORDER, LOWER_BORDER, EPSILON, SIGMA}; 
};
#endif // _WORLDLC_HPP
// vim:set et sts=4 ts=4 sw=4 ai ci cin:
