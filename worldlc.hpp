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

    /* @brief Overwrite existing World::readParameter to handle new files
     *
     * parameter file example TODO: rewrite code example
     * \code
     * \endcode
     *
     * @param filename filename of the parameter file
     */
    void read_Parameter(const std::string &filename);

    // Value-Defintions of the different String values
    // needed to be implemented again, because enum is not extandable
    enum Option { NAME, DELTA_T, T_END, LENGTH, UPPER_BORDER, LOWER_BORDER, EPSILON, SIGMA, CELLRCUT}; 
   
   /// cells
   std::vector<Cell> cells;
   /// Number of cells in every dimension
   int cell_N[DIM];
   /// length of cells
   real cell_length[DIM];
   /// r_cut used for calculation of the cell length
   real cell_r_cut;
};
#endif // _WORLDLC_HPP
// vim:set et sts=4 ts=4 sw=4 ai ci cin:
