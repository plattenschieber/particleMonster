#include "observerxyz.hpp"

ObserverXYZ::ObserverXYZ(WorldLC &_W) : Observer(_W), W(_W)
{
    std::ostringstream procID;
    procID << W.s.myrank;
    // open xyz file
    std::string coordFilename = "log/" + W.name + "_pid" + procID.str () + ".xyz";
    // open file, overwrite existing files, take no prisioners
    coordinates.open(coordFilename.c_str());
    if ( coordinates.is_open() )
        // and tell the world
        std::cout << "Opened " << coordFilename << " for writing." << std::endl;
}


ObserverXYZ::~ObserverXYZ()
{
    // close the coordinates file
    if ( coordinates.is_open() )
        coordinates.close();
    std::cout << "Everything closed properly" << std::endl;
}

void ObserverXYZ::outputCoordinates ()
{
    // write size and actual time of our world W
    xyz << W.nParticles << std::endl << "Time: " << W.t << std::endl;
    // get out every particle to satisfy the xyz format
    for (std::vector<Cell>::const_iterator i = W.cells.begin(); i < W.cells.end(); i++)
    {
        for (std::list<Particle>::const_iterator j = i->particles.begin(); j != i->particles.end(); j++)
        {
            // each particle should be an H-atom. At least now...
            xyz << "H\t";
            // particle j is located in a DIM-dimensional space
            for (unsigned int d=0; d<DIM; d++)
                // get it out, seperated with tabulars
                xyz << j->x[d] << "\t";
            // newline at end of each particle
            xyz << std::endl;
        }
    }
}

void ObserverXYZ::notify()
{
    // write statistics and coordinates
    Observer::notify();
    // write the xyz-format
    outputXYZ();
}



// vim:set et sts=4 ts=4 sw=4 ai ci cin:
