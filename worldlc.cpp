#include "worldlc.hpp"
#include "world.hpp"
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
    std::string line, option, tmp;
    
    // read file till eof
    while (parfile.good())
    {
        // read line from file
        getline(parfile,line);
        // create a string stream
        std::stringstream strstr;
        // put line into string stream
        strstr << line;
        // read option and value from stringstream
        strstr >> option;
        // push next read value, with internal converter of string stream, into the propper place 
        if (mapOptions[option] == CELLRCUT)
	    	strstr >> cell_r_cut;
    // close file
    parfile.close();
    // TODO: ADD Code to handle cell_r_cut and cell dimensions... 
}

void WorldLC::read_Particles(const std::string &filename)
{
    // call the base function
    World::read_Particles(filename);
    // TODO: Write every particle into it's cell
}

// vim:set et sts=4 ts=4 sw=4 ai ci cin: