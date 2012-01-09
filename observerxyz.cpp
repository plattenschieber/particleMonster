#include "observerxyz.hpp"

ObserverXYZ::ObserverXYZ(WorldLC &_W) : Observer(_W), W(_W)
{

    // open xyz file
    std::string xyzFilename = "log/" + W.name + ".xyz";
    // open file, overwrite existing files, take no prisioners
    xyz.open(xyzFilename.c_str());
    if ( xyz.is_open() )
        // and tell the world
        std::cout << "Opened " << xyzFilename << " for writing." << std::endl;
    
    // open coordinates file
    std::string coordFilename = "log/" + W.name + ".coordinates";
    // open file, overwrite existing files, take no prisioners
    coordinates.open(coordFilename.c_str());
    if ( coordinates.is_open() )
        // and tell the world
        std::cout << "Opened " << coordFilename << " for writing." << std::endl;
}


ObserverXYZ::~ObserverXYZ()
{
    // close the coordinates file
    if ( xyz.is_open() )
    {
        xyz.close();
        std::cout << "Closed xyz" << std::endl;
    }
    // close the coordinates file
    if ( coordinates.is_open() )
    {
        coordinates.close();
        std::cout << "Closed coordinates" << std::endl;
    }

    std::cout << "Everything closed properly" << std::endl;
}

void ObserverXYZ::outputCoordinates()
{
    // write updating time
    coordinates << ((WorldLC)W).t << "\t";
    // run over each particle,...
    for (std::vector<Cell>::const_iterator i = W.cells.begin(); i < W.cells.end(); i++)
        // each cell...
        for (std::list<Particle>::const_iterator j = i->particles.begin(); j != i->particles.end();	j++)
            // ...and each of it's dimensions
            for (unsigned int d=0; d<DIM; d++)
                // get it out, seperated with tabulars
                coordinates << j->x[d] << "\t";
    // end of line
    coordinates << std::endl;

}

void ObserverXYZ::outputXYZ()
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
