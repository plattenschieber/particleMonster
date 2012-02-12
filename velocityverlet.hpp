#ifndef _VELOCITYVERLET_HPP
#define _VELOCITYVERLET_HPP

#include "timediscretization.hpp"

/**
 * @brief Implementation of the Velocity Verlet Algorithm
 */
class VelocityVerlet : public TimeDiscretization {
public:
    /**
     * @brief constructor
     *
     * @param _W world configuration
     * @param _Pot potential used for force calculation
     * @param _O Observer of the simulation
     */
    VelocityVerlet(World& _W, Potential& _Pot, ParObserverXYZ &_O);
    
    /**
     * @brief constructor
     *
     * This is an example for Constructor overloading. If you have read until
     * here you can use the other constructor and change the blatt1 main function.
     *
     * @param _W world configuration
     * @param _Pot potential used for force calculation
     * @param _O Observer of the simulation
     */
    VelocityVerlet(World& _W, Potential* _Pot, ParObserverXYZ &_O);
    
    /**
     * @brief run a single timestep
     *
     * @param delta_t length of the timestep
     */
    void timestep(real delta_t);
    
    /**
     * @brief run the simulation
     */
    virtual void simulate();
    
    /**
     * @brief calculates the forces affecting the particles at the current time
     */
    virtual void compF();
    
    /**
     * @brief calculates the new velocity of the particles
     */
    virtual void updateV();
    
    /**
     * @brief calculate the new position of all particles according to their velocity
     */
    virtual void updateX();
    
    /**
     * @brief updates the average energy of the last 100 steps
     */
    void updateAverage();

    /**
     * @brief check every particle if it's still inside our word. Open borders will allow, energy loss in our system!
     */
    virtual void handleBorders();
    
protected:
    // data structures inherited from TimeDiscretization
    
private:
    VelocityVerlet();
};

#endif // _VELOCITYVERLET_HPP
// vim:set et sts=4 ts=4 sw=4 ai ci cin:
