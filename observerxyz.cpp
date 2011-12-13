#include "observerxyz.hpp"

ObserverXYZ::ObserverXYZ(WorldLC &_W) : Observer(_W), W(_W)
{

    // open xyz file
    std::string xyz_filename = W.name + ".xyz";
    // open file, overwrite existing files, take no prisioners
    xyz.open(xyz_filename.c_str());
    if ( xyz.is_open() )
        // and tell the world
        std::cout << "Opened " << xyz_filename << " for writing." << std::endl;
    
    // open coordinates file
    std::string coordinates_filename = W.name + ".coordinates";
    // open file, overwrite existing files, take no prisioners
    coordinates.open(coordinates_filename.c_str());
    if ( coordinates.is_open() )
        // and tell the world
        std::cout << "Opened " << coordinates_filename << " for writing." << std::endl;
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

void ObserverXYZ::output_coordinates()
{
    // write updating time
    coordinates << ((WorldLC)W).t << "\t";
    // run over each particle...
    for (std::vector<Cell>::const_iterator i = ((WorldLC)W).cells.begin(); i < ((WorldLC)W).cells.end(); i++)
	// TODO: Change in other place, too
	for (std::vector<Particle>::const_iterator j = i->particles.begin(); j < i->particles.begin();	j++)
	{
 	    // ...and each of it's dimensions
            for (unsigned int d=0; d<DIM; d++)
                // get it out, seperated with tabulars
                coordinates << j->x[d] << "\t";
            // end of line
	    coordinates << std::endl;
	}
}

void ObserverXYZ::output_xyz()
{
    // write size and actual time of our world W
    xyz << W.nParticles << std::endl << "Time: " << W.t << std::endl;
    // get out every particle to satisfy the xyz format
    for (std::vector<Cell>::const_iterator i = W.cells.begin(); i < W.cells.end(); i++)
    {
        for (std::vector<Particle>::const_iterator j = i->particles.begin(); j < i->particles.end(); j++)    
    	{
            // each particle should be an H-atom. At least now... 
            xyz << "H\t";
            std::cout <<  "H\t";
            // particle j is located in a DIM-dimensional space
            for (unsigned int d=0; d<DIM; d++)
            {
                // get it out, seperated with tabulars
                xyz << j->x[d] << "\t";
                std::cout << j->x[d] << "\t";
            }
            xyz << std::endl;
            std::cout << std::endl;
	    }
	    // new line at end of one timestep
        if(i == W.cells.end()-1); 
            //xyz << std::endl;
    }
}

void ObserverXYZ::notify()
{
    // write statistics and coordinates 
    Observer::notify();
    // write the xyz-format
   // std::cout << "NOTIFYXYZ" << std::endl;
    output_xyz();
}



// vim:set et sts=4 ts=4 sw=4 ai ci cin:
