#ifndef _WORLDLC_HPP
#define _WORLDLC_HPP

#include "world.hpp"
#include "cell.hpp"
#include "defines.hpp"
#include "particle.hpp"
#include <mpi.h>
#include <vector>
#include <list>
#include <string>
#include <iostream>
#include <fstream>
#include <map>
/**
 * @brief the world class holds all information of the simulation environment
 */
class WorldLC : public World {
public:
    // Standard Ctor
    WorldLC();

    /* @brief Overwrite existing World::readParameter to handle new files
     *
     * parameter file example TODO: rewrite code example
     * \code
     * \endcode
     *
     * @param filename filename of the parameter file
     */
    virtual void readParameter (const std::string &filename);
    // TODO: Add some comment here
    virtual void readParticles (const std::string &filename);

    // get the global cell number of a given particle
    virtual int getCellNumber (const Particle &p);

    // Value-Defintions of the different String values
    // needed to be implemented again, because enum is not extandable
    enum Option { NAME, DELTA_T, T_END, LENGTH, UPPERBORDER, LOWERBORDER, EPSILON, SIGMA, CELLRCUT};
    // Map to associate the strings with the enum values
    std::map<std::string, WorldLC::Option> mapOptions;
    /// a vector of cells with particle lists
    std::vector<Cell> cells;
    /// Number of cells in every dimension
    int nCells[DIM];
    /// length of cells
    real cellLength[DIM];
    /// r_cut used for calculation of the cell length
    real cell_r_cut;
};

/**
 * @brief a ostream operator for the WorldLC class
 *
 * @param os stream object
 * @param W the world
 *
 * @return resulting stream object
 */
std::ostream& operator << (std::ostream& os, WorldLC& W);


#endif // _WORLDLC_HPP
// vim:set et sts=4 ts=4 sw=4 ai ci cin:
