#include "observer.hpp"

Observer::Observer(World &_W) : W(_W)
{

    // open xyz file
    std::string xyz_filename = "/xyz/" + W.name + ".xyz";
    // open file, overwrite existing files, take no prisioners
    xyz.open(xyz_filename.c_str());
    if ( xyz.is_open() )
        // and tell the world
        std::cout << "Opened " << xyz_filename << " for writing." << std::endl;
    
    // open statistics file
    std::string statistics_filename = "/statistics/" + W.name + ".statistics";
    // open file, overwrite existing files, take no prisioners
    statistics.open(statistics_filename.c_str());
    if ( statistics.is_open() )
        // and tell the world
        std::cout << "Opened " << statistics_filename << " for writing." << std::endl;
    
    // open coordinates file
    std::string coordinates_filename = "/coordinates/" + W.name + ".coordinates";
    // open file, overwrite existing files, take no prisioners
    coordinates.open(coordinates_filename.c_str());
    if ( coordinates.is_open() )
        // and tell the world
        std::cout << "Opened " << coordinates_filename << " for writing." << std::endl;
}


Observer::~Observer()
{
    // close the coordinates file
    if ( xyz.is_open() )
    {
        xyz.close();
        std::cout << "Closed xyz" << std::endl;
    }
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

void Observer::output_statistics()
{
    // write statistics into the filestream, seperated with tabulars
    statistics
        << W.t << "\t" 
        << W.e_pot << "\t"
        << W.e_kin << "\t"
        << W.e_kin + W.e_pot // observe conservation of energy 
        << std::endl;
}

void Observer::output_coordinates()
{
    // write updating time
    coordinates << W.t << "\t";
    // run over each particle...
    for (std::vector<Particle>::const_iterator i = W.particles.begin(); i != W.particles.end(); ++i)
        // ...and each of it's dimensions
        for (unsigned int d=0; d<DIM; d++)
            // get it out, seperated with tabulars
            coordinates << i->x[d] << "\t";
    // end of line
    coordinates << std::endl;
}

void Observer::output_xyz()
{
    // write size and actual time of our world W
    xyz << W.particles.size() << std::endl << "Time: " << W.t << std::endl;
    // get out every particle to satisfy the xyz format
    for (std::vector<Particle>::const_iterator i = W.particles.begin(); i != W.particles.end(); ++i)
    {
        // each particle should be an H-atom. At least now... 
        xyz << "H\t";
        // particle i is located in a DIM-dimensional space
        for (unsigned int d=0; d<DIM; d++)
            // get it out, seperated with tabulars
            xyz << i->x[d] << "\t";
        // new line only after the last particle
        if(i != W.particles.end()-1) xyz << std::endl;
    }
    // new line at end of one timestep
    xyz << std::endl;
    
}

void Observer::notify()
{
    // write statistics 
    output_statistics();
    // write the coordinates of our particles
    output_coordinates();
    // write the xyz-format
    output_xyz();
}



// vim:set et sts=4 ts=4 sw=4 ai ci cin:
