#include "worldlc.hpp"
#include "world.hpp"
#include "cell.hpp"
#include <stdexcept>
#include <sstream>
#include <map>


// ctor, which calls World::World()   
WorldLC::WorldLC() {
    // we do need another mapOption
    mapOptions["cell_r_cut"] = CELLRCUT;
}


void WorldLC::read_Parameter(const std::string &filename)
{
    // call the base function 
    World::read_Parameter(filename); 
    // create input filestream
    std::ifstream parfile(filename.c_str());
    // check if file is open
    if (!parfile.is_open())
        throw std::runtime_error("read_Parameter(): Can't open file '" + filename + "' for reading.");
    
    // helper strings
    std::string line, option;
    
    // read file till eof
    while (parfile.good())
    {
        // read line from file
        getline(parfile,line);
        // create a string stream
        std::stringstream strstr;
        // put line into string stream
        strstr << line;
        // read option from stringstream
        strstr >> option;
        // TODO: Add eps and sigma!
        // push next read value, with internal converter of string stream, into the propper place 
        if (mapOptions[option] == CELLRCUT)
	    	strstr >> cell_r_cut;
    }
    // close file
    parfile.close();
    
    // #of all cells
    int nCells = 1;

    for (int d=0; d<DIM; d++) 
    {
        std::cout << std::endl << "length: " << length[d] << "cell r cut: " << cell_r_cut << std::endl;
        // #cells in dimension = floor(length per cell-cutlength)
        cell_N[d] = (int)(length[d]/cell_r_cut);
        // cell length = world length per #cells
        cell_length[d] = length[d]/cell_N[d];

        nCells *= cell_N[d];
        // DEBUG
        std::cout << "#Cell: " << cell_N[d] << "\tCelllength: " << cell_length[d] << std::endl;
    }
    // DEBUG 
    std::cout << "#Cells: " << nCells << std::endl;

    // hold on computer and give us all what you have! 
    std::cin.clear(); 
    std::cin.ignore(INT_MAX, '\n');  
    std::cout << "Press Return..." << std::endl;  
    std::cin.get();


    
    
    // insert empty cells
    for (int i=0; i<nCells; i++);
        cells.push_back(Cell());
}

void WorldLC::read_Particles(const std::string &filename)
{
    // call the base function
    World::read_Particles(filename);
    // Write every particle into it's belonging cell
    for (std::vector<Particle>::iterator i = particles.begin(); i < particles.end(); i++)
    {
        // add particle to right cell...
        // getCellNumber(i) gives belonging cellnumber, push into this cell our actual particle i: particles[i-particles.begin()]
        cells[getCellNumber(i)].particles.push_back(particles[i-particles.begin()]);
            //particles.push_back(i);
        // ...and remove it from our particle list
        particles.erase(i);

    }
}
int WorldLC::getCellNumber(const std::vector<Particle>::iterator i) 
{
    int tmp[3] = {0,0,0};
    for (int d=0; d<DIM; d++)
    {
        tmp[d] = i->x[d] * cell_N[d] / length[d];
        // DEBUG
        std::cout << tmp[d] << "\tx[" << d << "]: " << i->x[d] << "\t "<< cell_N[d] / length[d] << std::endl;
    }
    std::cout << "Corresponding Index: " << J(tmp,cell_N) << std::endl;
    return J(tmp,cell_N);

}
// vim:set et sts=4 ts=4 sw=4 ai ci cin:
