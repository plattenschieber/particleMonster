#include <cstdlib>
#include <iostream>
#include <cstdio>
#include <cmath>
#include <mpi.h>


#include "world.hpp"
#include "worldlc.hpp"
#include "gravitypotential.hpp"
#include "ljpotential.hpp"
#include "velocityverlet.hpp"
#include "velocityverletlc.hpp"
#include "observerxyz.hpp"

// const char not available for mpi
int main(int argc, char *argv[]) {

    int myrank, numprocs;

    // Initialize MPI.
     MPI::Init (argc, argv);

     // Get the number of processes.
     numprocs = MPI::COMM_WORLD.Get_size();
     // Get the individual process ID.
     myrank = MPI::COMM_WORLD.Get_rank();

    // check arguments
    if (argc < 2) {
        std::cerr << "error: missing arguments" << std::endl;
        std::cerr << "usage: " << std::endl
                  << "\t" << argv[0] << " parameterfile particledatafile" << std::endl;
        return EXIT_FAILURE;
    }

    // instanciate Potentials
    LJPotential LPot;

    // create World
    WorldLC W;
    // debub propose, show 2nd argument
    std::cout << argv[2] << std::endl;

    // read Parameters
    W.readParameter(argv[1]);

    // read Particles
    W.readParticles(argv[2]);

    // print World configuration after Building it up
    std::cout << W << std::endl;

    // create the Observer
    ObserverXYZ O(W);

    // instanciate timediscretization
    // remark: & is used to get the address of Pot // reremark: removed the &
    VelocityVerletLC verletLC(W, LPot, O);
    std::cout << "WHAT THE FUCK" << std::endl;

    // run the simulation
    verletLC.simulate();

    return EXIT_SUCCESS;
}

// vim:set et sts=2 ts=2 sw=2 ai ci cin:
