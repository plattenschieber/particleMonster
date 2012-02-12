#include "worldlc.hpp"
#include "world.hpp"
#include "cell.hpp"
#include <stdexcept>
#include <sstream>
#include <map>
#include <cmath>

// ctor, which calls World::World()   
WorldLC::WorldLC() : cell_r_cut(2.5) {
    // we do need another mapOption
    mapOptions["cell_r_cut"] = CELLRCUT;
    mapOptions["num_procs"] = NUMPROCS;
}


void WorldLC::readParameter(const std::string &filename)
{
    // call the base function
    World::readParameter(filename);
    // create input filestream
    std::ifstream parfile(filename.c_str());
    // check if file is open
    if (!parfile.is_open())
        throw std::runtime_error("readParameter(): Can't open file '" + filename + "' for reading.");
    
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
        // push next read value, with internal converter of string stream, into the propper place
        if (mapOptions[option] == CELLRCUT)
            strstr >> cell_r_cut;
        if (mapOptions[option] == NUMPROCS)
            for (int d=0; d<DIM; d++)
                strstr >> nProcs[d];
    }
    // close file
    parfile.close();

    // Calc #cells
    for (int d=0; d<DIM; d++)
    {
        // #cells in dimension = floor(length per cell-cutlength)
        nCells[d] = (int)(worldLength[d]/cell_r_cut);
        // cell length = world length per #cells
        cellLength[d] = worldLength[d]/nCells[d];
    }

    // calc #of all cells incl. bordures
    int numCells = 1;
    for (int d=0; d<DIM; d++)
        numCells *= cellLength[d];
    // insert empty cells
    //for (int i=0; i<numCells; i++)
    //    cells.push_back(Cell());
    // resize vector size
    cells.resize (numCells);

    std::cout << "END OF readParameter()" << std::endl << this;
}


void WorldLC::readParticles(const std::string &filename)
{
    // call the base function
    World::readParticles(filename);
    nParticles = particles.size();
    // Write every particle into it's belonging cell
    for (std::list<Particle>::iterator i = particles.begin(); i != particles.end(); i++)
    {
        std::cout << "Push and erase particles[" << i->ID << "] " << std::endl;
        // add particle to right cell...
        // getCellNumber(i) gives belonging cellnumber, push into this cell our actual particle i
        cells[getCellNumber(*i)].particles.push_back(*i);
    }
    std::cout << "***************************************************" << std::endl
              << "FINISHED READING PARTICLES - NOW CLEAR PARTICLES" << std::endl
              << "*************************************************** \n" << std::endl;
    // Clear all particles in our World Vector
    particles.clear();
}


int WorldLC::getCellNumber(const Particle &p)
{
    // temporary array
    int tmp[3] = {0,0,0};
    for (int d=0; d<DIM; d++)
    {
        // handle particle outside the world failure
        if (p.x[d] < 0 || p.x[d] > worldLength[d])
        {
            std::cerr << "<------- FAILURE ------->" << std::endl;
            std::cerr << "Particle left faulty the lower/upper border - in getCellNumber() -> x["
                      << d << "] = " << p.x[d] << std::endl;
            exit(EXIT_FAILURE);
        }
        // compute global cell number
        tmp[d] = (int) floor(p.x[d] / cellLength[d]) % nCells[d];
    }
    // return corresponding cell Number
    return J(tmp,nCells);
}


std::ostream& operator << (std::ostream& os, WorldLC& W) 
{
    // Get out some information about the world and it's SubDomain
    os << W.name << " Dim: " << DIM << " t: " << W.t << " delta_t: " << W.delta_t << " t_end: " << W.t_end
       << " cell_r_cut: " << W.cell_r_cut << std::endl;
    os << "Length: ";
    for (int d=0; d<DIM; d++)
        os << W.worldLength[d] << " ";
    os << std::endl << "Upper borders: ";
    for (int d=0; d<DIM; d++)
        os << W.upper_border[d] << " ";

    os << std::endl << "Lower borders: ";
    for (int d=0; d<DIM; d++)
        os << W.lower_border[d] << " ";
    os << std::endl << "Border types: 0 - leaving, 1 - periodic, 2 - unknown" << std::endl;

    os << std::endl << "Present particles in all " << W.cells.size() << " Cells." << std::endl;

    // roll over each Cell
    for (std::vector<Cell>::iterator i = W.cells.begin(); i < W.cells.end(); i++)
    {
        // if there are some particles in this cell, show it off
        if (i->particles.size() > 0)
        {
            // we are now in cell# x
            os << "cell[" << i-W.cells.begin() << "] = " ;
            // get out every particle in this cell
            for (std::list<Particle>::iterator j = i->particles.begin(); j != i->particles.end(); j++)
            {
                // Particle Number
                os <<  j->ID << ": ";
                //  get out the particles location
                for (unsigned int d=0; d<DIM; d++)
                    os << j->x[d] << "\t";
                // tabulator after each particle
                os << "\t";
            }
            // newline after each cell
            os << std::endl;
        }
    }
    return os;
}

// vim:set et sts=4 ts=4 sw=4 ai ci cin:
