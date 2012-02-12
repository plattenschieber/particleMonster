#include <cstdlib>
#include <iostream>
#include <cstdio>
#include <cmath>
#include <mpi.h>


#include "parobserverxyz.hpp"
#include "world.hpp"
#include "worldlc.hpp"
#include "gravitypotential.hpp"
#include "ljpotential.hpp"
#include "velocityverlet.hpp"
#include "velocityverletlc.hpp"
#include "observerxyz.hpp"
#include "subdomain.hpp"

int main(int argc, char *argv[]) {
    // Initialize MPI.
     MPI::Init (argc, argv);

    // check arguments
    if (argc != 3) {
        if(argc < 3)
        {
            std::cerr << "error: missing arguments" << std::endl;
            MPI::Finalize();
            return EXIT_FAILURE;
        }
        else if (argc > 3)
        {
            std::cerr << "error: too many arguments" << std::endl
                      <<  "usage: " << std::endl
                      << "\t" << argv[0] << " parameterfile particledatafile" << std::endl;
            MPI::Finalize();
            return EXIT_FAILURE;
        }
    }

    // instanciate Potentials
    LJPotential LPot;

    // create World
    SubDomain S;
    // debub propose, show 2nd argument
    std::cout << argv[2] << std::endl;

    // read Parameters
    S.readParameter(argv[1]);

    // read Particles
    S.readParticles(argv[2]);

    WorldLC &W = S;
    // print World configuration after Building it up
    std::cout << S << std::endl;
    std::cout << W << std::endl;

    // create the Observer
    ParObserverXYZ O(S);

    // instanciate timediscretization
    // remark: & is used to get the address of Pot // reremark: removed the &
    VelocityVerletLC verletLC(S, LPot, O);

    // run the simulation
    verletLC.simulate();

    // terminate MPI
    MPI::Finalize();

    return EXIT_SUCCESS;
}

// vim:set et sts=2 ts=2 sw=2 ai ci cin:
