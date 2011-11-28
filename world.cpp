#include "world.hpp"
#include <stdexcept>
#include <sstream>
#include <map>


// intializing via intializer list - order is set by definings in .hpp; 
World::World() : name("unknown"),t(0),delta_t(0.0),t_end(0.0),e_kin(0.0),e_pot(0.0),e_tot(0.0)
{    
    //setting map mapOptions for switchcase
    mapOptions["name"] = NAME;
    mapOptions["delta_t"] = DELTA_T;  
    mapOptions["t_end"] = T_END;
    mapOptions["length"] = LENGTH;
    mapOptions["upper_border"] = UPPER_BORDER;
    mapOptions["lower_border"] = LOWER_BORDER;
    mapOptions["sigma"] = SIGMA;
    mapOptions["epsilon"] = EPSILON;
}

virtual void World::read_Parameter(const std::string &filename)
{
  
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
        switch(mapOptions[option])
        {
            case DELTA_T:
                strstr >> delta_t;
                break;
            case T_END:
                strstr >> t_end;
                break;
            case NAME:
                strstr >> name;
                break;
            case LENGTH:
                for (int i=0; i<DIM; i++) {
                    strstr >> length[i];
                } 
                break;
            case UPPER_BORDER:
                strstr >> tmp;
                for (int i=0; i<DIM; i++)
                {
                    if (tmp == "leaving") upper_border[i] = leaving;   
		    else if (tmp == "periodic") upper_border[i] = periodic;
                    else upper_border[i] = unknown;
                }
                break;
            case LOWER_BORDER:
                strstr >> tmp;
                for (int i=0; i<DIM; i++)
                {
                    if (tmp == "leaving") lower_border[i] = leaving;   
		    else if (tmp == "periodic") upper_border[i] = periodic;
                    else lower_border[i] = unknown;
                }
                break;
        }
    }
    // close file
    parfile.close();
}

virtual void World::read_Particles(const std::string &filename)
{
    // create input filestream
    std::ifstream parfile(filename.c_str());
    // check if file is open
    if (!parfile.is_open())
        throw std::runtime_error("read_Particles(): Can't open file '" + filename + "' for reading.");
    
    // helper strings
    std::string line;
    Particle tmpparticle;
    
    // read file till eof
    while (parfile.good())
    {
        // read line from file
        getline(parfile,line);
        if (line=="")
            break;
        // create a string stream
        std::stringstream strstr;
        // push line into string stream
        strstr << line;    
        
        // grab ID
        strstr >> tmpparticle.ID;        
        // save the particles mass
        strstr >> tmpparticle.m; 
        // push dim read values into tmpparticles position
        for(int i=0; i<DIM; i++)
            strstr >> tmpparticle.x[i];
        // push dim read values into tmpparticles velocity
        for(int i=0; i<DIM; i++)
            strstr >> tmpparticle.v[i];        
        // assure integrity of the force of our tmpparticle
        for(int i=0; i<DIM; i++)
            tmpparticle.F[i] = tmpparticle.F_old[i] = 0.0;
        
        // add the new particle to our particles
        particles.push_back(tmpparticle);
        // resets tmpparticle for next 
        tmpparticle.clear();
    }
    // close file
    parfile.close();
}

std::ostream& operator << (std::ostream& os, World& W) {
    os << "t=" << W.t << " delta_t=" << W.delta_t << " t_end=" << W.t_end
    << " Number of Particles=" << W.particles.size();
    for (unsigned int i=0; i<W.particles.size(); i++) {
        os << "particle[" << i << "] = (";
        for (unsigned int d=0; d<DIM; d++)
            os << W.particles[i].x[d] << ", ";
        os << ")";
    }
    return os;
}
// vim:set et sts=4 ts=4 sw=4 ai ci cin:
