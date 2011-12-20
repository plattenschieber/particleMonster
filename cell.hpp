#ifndef _CELL_HPP
#define _CELL_HPP

#include "defines.hpp"
#include <list>
#include "particle.hpp"

/**
 * @brief Particle data
 *
 * this class contains the particle data
 */
class Cell {
public:
    /// The cell contains particles ...
    std::list<Particle> particles;
};

#endif // _CELL_HPP
// vim:set et sts=4 ts=4 sw=4 ai ci cin:
