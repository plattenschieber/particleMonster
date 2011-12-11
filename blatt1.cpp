#include <cstdlib>
#include <iostream>
#include <cstdio>
#include <math.h>

#include "world.hpp"
#include "gravitypotential.hpp"
#include "ljpotential.hpp"
#include "velocityverlet.hpp"
#include "observer.hpp"

int main2(int argc, const char *argv[]) {
    
    // check arguments
    if (argc < 2) {
        std::cerr << "error: missing arguments" << std::endl;
        std::cerr << "usage: " << std::endl
        << "\t" << argv[0] << " parameterfile particledatafile" << std::endl;
        return EXIT_FAILURE;
    }
    
    // instanciate Potentials
    GravityPotential Pot;
    LJPotential LPot;
    
    // create World
    World W;
    // debub propose, show 2nd argument
    std::cout << argv[2] << std::endl;
    
    // read Parameters
    W.read_Parameter(argv[1]);
    
    // read Particles
    W.read_Particles(argv[2]);
    
    // print World configuration after Building it up
    std::cout << W << std::endl;
    
    // create the Observer
    Observer O(W);
    
    // instanciate timediscretization 
    // remark: & is used to get the address of Pot // reremark: removed the & 
    VelocityVerlet Verlet(W, LPot, O);
    
    // run the simulation
    Verlet.simulate();
    
    return EXIT_SUCCESS;
}

// vim:set et sts=2 ts=2 sw=2 ai ci cin:
