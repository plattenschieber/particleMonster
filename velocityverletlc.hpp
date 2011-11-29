#ifndef _VELOCITYVERLETLC_HPP
#define _VELOCITYVERLETLC_HPP

#include "worldlc.hpp"
#include "velocityverlet.hpp"

/**
 * @brief Implementation of the Velocity Verlet Algorithm by means of our cell structure
 * Inheritance from VelocityVerlet, because we don't need to reimplement constructor, simulate, ... 
 */
class VelocityVerletLC : public VelocityVerlet {
public:
    /**
     * @brief constructor
     *
     * @param _W world configuration
     * @param _Pot potential used for force calculation
     * @param _O Observer of the simulation
     */
    VelocityVerletLC(WorldLC& _W, Potential& _Pot, Observer &_O);
    
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
    VelocityVerletLC(WorldLC& _W, Potential* _Pot, Observer &_O);
    
    WorldLC W;
    
    /**
     * @brief calculates the forces affecting the particles at the current time
     */
    void comp_F();
    
    /**
     * @brief calculates the new velocity of the particles
     */
    void update_V();
    
    /**
     * @brief calculate the new position of all particles according to their velocity
     */
    void update_X();
    
    /**
     * @brief check every particle if it's still inside our word. Open borders will allow, energy loss in our system!
     */
    void handle_borders();
    
protected:
    // data structures inherited from TimeDiscretization
    
private:
    VelocityVerletLC();
};

#endif // _VELOCITYVERLETLC_HPP
// vim:set et sts=4 ts=4 sw=4 ai ci cin:
