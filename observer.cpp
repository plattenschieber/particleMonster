#include "observer.hpp"

Observer::Observer(World &_W) : W(_W)
{
    // open statistics file
    std::string statistics_filename = "log/" + W.name + ".statistics";
    // open file, overwrite existing files, take no prisioners
    statistics.open(statistics_filename.c_str());
    if ( statistics.is_open() )
        // and tell the world
        std::cout << "Opened " << statistics_filename << " for writing." << std::endl;
    
    // open coordinates file
    std::string coordinates_filename = "log/" + W.name + ".coordinates";
    // open file, overwrite existing files, take no prisioners
    coordinates.open(coordinates_filename.c_str());
    if ( coordinates.is_open() )
        // and tell the world
        std::cout << "Opened " << coordinates_filename << " for writing." << std::endl;
}


Observer::~Observer()
{
    // close the statistics file
    if ( statistics.is_open() )
    {
        statistics.close();
        std::cout << "Closed statistics" << std::endl;
    }
    // close the coordinates file
    if ( coordinates.is_open() )
    {
        coordinates.close();
        std::cout << "Closed coordinates" << std::endl;
    }

    std::cout << "Everything closed properly" << std::endl;
}

void Observer::outputStatistics()
{
    // write statistics into the filestream, seperated with tabulars
    statistics
            << W.t << "\t"
            << W.e_pot << "\t"
            << W.e_kin << "\t"
            << W.e_kin + W.e_pot << "\t" // observe conservation of energy
            << W.e_avg // observe energy conservation of last 100 steps
            << std::endl;
}

void Observer::outputCoordinates()
{
    // write updating time
    coordinates << W.t << "\t";
    // run over each particle (use const where you can)...
    for (std::list<Particle>::const_iterator i = W.particles.begin(); i != W.particles.end(); ++i)
        // ...and each of it's dimensions
        for (unsigned int d=0; d<DIM; d++)
            // get it out, seperated with tabulars
            coordinates << i->x[d] << "\t";
    // end of line
    coordinates << std::endl;
}

void Observer::notify()
{
    // write statistics
    outputStatistics();
    // write the coordinates of our particles
    outputCoordinates();
}



// vim:set et sts=4 ts=4 sw=4 ai ci cin:
