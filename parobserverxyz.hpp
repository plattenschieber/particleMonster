#ifndef _PAROBSERVERXYZ_HPP
#define _PAROBSERVERXYZ_HPP

#include "subdomain.hpp"
#include "observerxyz.hpp"
#include <iostream>
#include <fstream>

/**
 * @brief This Observer will controll the coordinates filenames and limit statistics output to first proc
 */
class ParObserverXYZ : public ObserverXYZ  {
public:
    /**
     * @brief constructor
     *
     * opens and creates the files written during observation
     *
     * @param _W
     */
    ParObserverXYZ(SubDomain& _S);

    /**
     * @brief destructor
     *
     * closes the files written during the obervation
     */
    ~ParObserverXYZ();

    /*
    * @brief notify the observer that the world has changed
    */
    virtual void notify ();


    
protected:
    /// The world we are observing
    SubDomain &S;

private:
    /// Disabled Standard Constructor
    ParObserverXYZ();
};

#endif // _PAROBSERVERXYZ_HPP
// vim:set et sts=4 ts=4 sw=4 ai ci cin:
