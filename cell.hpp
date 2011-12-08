#ifndef _CELL_HPP
#define _CELL_HPP

#include "defines.hpp"
#include <vector>

/**
 * @brief Particle data
 *
 * this class contains the particle data
 */
class Cell {
public:
    // standard ctor 
    Cell();
    /// The cell contains particles ...
    std::vector<Particle> particles;
    // DEBUG: 
    //int position; //int groesse
};

#endif // _CELL_HPP
// vim:set et sts=4 ts=4 sw=4 ai ci cin:
