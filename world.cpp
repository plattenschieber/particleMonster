#include "world.hpp"
#include <stdexcept>
#include <sstream>
#include <map>


    // intializing via intializer list - order is set by definings in .hpp; 
    World::World() : name("unknown"),t(0),deltaT(0.0),tEnd(0.0),eKin(0.0),ePot(0.0),eTot(0.0),sigma(1.0),epsilon(1.0)
    {    
        //setting map mapOptions for switchcase
        mapOptions["name"] = NAME;
        mapOptions["delta_T"] = DELTA_T;
        mapOptions["t_End"] = T_END;
        mapOptions["length"] = LENGTH;
        mapOptions["upper_border"] = UPPERBORDER;
        mapOptions["lower_border"] = LOWERBORDER;
        mapOptions["sigma"] = SIGMA;
        mapOptions["epsilon"] = EPSILON;
    }

    void World::readParameter(const std::string &filename)
    {
      
        // create input filestream
        std::ifstream parfile(filename.c_str());
        // check if file is open
        if (!parfile.is_open())
            throw std::runtime_error("readParameter(): Can't open file '" + filename + "' for reading.");
        
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
                    strstr >> deltaT;
                    break;
                case T_END:
                    strstr >> tEnd;
                    break;
                case NAME:
                    strstr >> name;
                    break;
                case EPSILON:
                    strstr >> epsilon;
                case SIGMA:
                    strstr >> sigma;
                case LENGTH:
                    for (int i=0; i<DIM; i++) {
                        strstr >> length[i];
                    } 
                    break;
                case UPPERBORDER:
                    strstr >> tmp;
                    for (int i=0; i<DIM; i++)
                    {
                        if (tmp == "leaving") upperBorder[i] = leaving;
                        else if (tmp == "periodic") upperBorder[i] = periodic;
                        else upperBorder[i] = unknown;
                    }
                    break;
                case LOWERBORDER:
                    strstr >> tmp;
                    for (int i=0; i<DIM; i++)
                    {
                        if (tmp == "leaving") lowerBorder[i] = leaving;
                        else if (tmp == "periodic") lowerBorder[i] = periodic;
                        else lowerBorder[i] = unknown;
                    }
                    break;
                // handle unknown options
                default:
                    std::cout << "'" << option << "' is an invalid option." << std::endl;
                    break;
             }
        }
        // close file
        parfile.close();
    }

    void World::readParticles(const std::string &filename)
    {
        // create input filestream
        std::ifstream parfile(filename.c_str());
        // check if file is open
        if (!parfile.is_open())
            throw std::runtime_error("readParticles(): Can't open file '" + filename + "' for reading.");
        
        // helper strings
        std::string line;
        Particle tmpparticle;
        
        // read file till eof
        while (parfile.good())
        {
            // read line from file
            getline(parfile,line);
            // break if there is a blank line, e.g. a newline at the eof 
            if (line=="")
                break;
            // create a string stream
            std::stringstream strstr;
            // push line into string stream
            strstr << line;    
            // grab ID:
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
            
            // add the new particle to our worlds' particles
            particles.push_back(tmpparticle);
            // resets tmpparticle for next 
            tmpparticle.clear();
    }
    // close file
    parfile.close();
}

std::ostream& operator << (std::ostream& os, World& W) {
    os << "t=" << W.t << " deltaT=" << W.deltaT << " tEnd=" << W.tEnd
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
