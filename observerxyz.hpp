#ifndef _OBSERVERXYZ_HPP
#define _OBSERVERXYZ_HPP

#include "worldlc.hpp"
#include "observer.hpp"
#include <iostream>
#include <fstream>

/**
 * @brief This Observer will look for our XYZ Data
 */
class ObserverXYZ : public Observer  {
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
    virtual ~ObserverXYZ();

    /*
    * @brief notify the observer that the world has changed
    */
    virtual void notify ();


    /**
     * @brief output coordinates of the particles in XYZ Format
     */
    virtual void outputCoordinates();
    
protected:
    /// The world we are observing
    WorldLC &W;
    /// coordiantes filestream
    std::ofstream coordinates;

private:
    /// Disabled Standard Constructor
    ObserverXYZ();
};

#endif // _OBSERVERXYZ_HPP
// vim:set et sts=4 ts=4 sw=4 ai ci cin:
