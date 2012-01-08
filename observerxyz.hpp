#ifndef _OBSERVERXYZ_HPP
#define _OBSERVERXYZ_HPP

#include "worldlc.hpp"
#include "observer.hpp"
#include <iostream>
#include <fstream>

/**
 * @brief This Observer will look for our XYZ Data
 */
class ObserverXYZ : Observer  {
public:
    /**
     * @brief constructor
     *
     * opens and creates the files written during observation
     *
     * @param _W
     */
    ObserverXYZ(WorldLC& _W);

    /**
     * @brief destructor
     *
     * closes the files written during the obervation
     */
    ~ObserverXYZ();
    
    /**
     * @brief notify the observer that the world has changed
     */
    void notify(); 
    
    /**
     * @brief generate a pymol readable file
     */
    void output_xyz();

    /** 
     * @brief output coordinates of the particles
     */
    void output_coordinates();
    
  

protected:
    /// The world we are observing
    WorldLC &W;
    /// filestream for xyz data
    std::ofstream xyz;
    /// coordiantes filestream
    std::ofstream coordinates;

private:
    /// Disabled Standard Constructor
    ObserverXYZ();
};

#endif // _OBSERVERXYZ_HPP
// vim:set et sts=4 ts=4 sw=4 ai ci cin:
