#include "parobserverxyz.hpp"

ParObserverXYZ::ParObserverXYZ(SubDomain &_S) : ObserverXYZ(_S), S(_S)
{
    std::ostringstream procID;
    procID << S.myrank;
    // open xyz file
    std::string coordFilename = "log/" + S.name + "_pid" + procID.str () + ".xyz";
    // open file, overwrite existing files, take no prisioners
    coordinates.open(coordFilename.c_str());
    if ( coordinates.is_open() )
        // and tell the world
        std::cout << "Opened " << coordFilename << " for writing." << std::endl;
}

void ParObserverXYZ::notify()
{
    // call output functions
    if(S.myrank==0)
        outputStatistics();
    outputCoordinates ();
}
ParObserverXYZ::~ParObserverXYZ()
{
    // close the coordinates file
    if ( coordinates.is_open() )
        coordinates.close();
    std::cout << "Everything closed properly in ParObserverXYZ" << std::endl;
}
// vim:set et sts=4 ts=4 sw=4 ai ci cin:
