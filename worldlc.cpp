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
        // TODO: Add eps and sigma!
        // push next read value, with internal converter of string stream, into the propper place
        if (mapOptions[option] == CELLRCUT)
            strstr >> cell_r_cut;
    }
    // close file
    parfile.close();
    

    // #of all cells
    int nCells = 1;
    // Calc #cells
    for (int d=0; d<DIM; d++)
    {
        // #cells in dimension = floor(length per cell-cutlength)
        cell_N[d] = (int)(length[d]/cell_r_cut);
        // cell length = world length per #cells
        cell_length[d] = length[d]/cell_N[d];

        nCells *= cell_N[d];
        // DEBUG
        //        std::cout << "World.length[" << d << "]=" << length[d] << "\tcell_r_cut=" << cell_r_cut
        //                  << "\t#Cells=" << cell_N[d] << "\tCelllength=" << cell_length[d] << std::endl;
    }
    
    // insert empty cells
    for (int i=0; i<nCells; i++)
        cells.push_back(Cell());

    // DEBUG
    //    std::cout << "#Cells: " << nCells << "\t" << cells.size() <<std::endl << std::endl;
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
        // getCellNumber(i) gives belonging cellnumber, push into this cell our actual particle i: particles[i-particles.begin()]

        cells[getCellNumber(i)].particles.push_back(*i);
    }
    std::cout << "***************************************************" << std::endl
              << "FINISHED READING PARTICLES - NOW CLEAR PARTICLES" << std::endl
              << "*************************************************** \n" << std::endl;
    // Clear all particles in our World Vector
    particles.clear();
}

int WorldLC::getCellNumber(const std::list<Particle>::iterator i)
{
    int tmp[3] = {0,0,0};
    //    // DEBUG Table
    //    std::cout << "Cell coordinate: " ;
    for (int d=0; d<DIM; d++)
    {
        if (i->x[d] < 0)
            return -1;
        tmp[d] = (int) floor(i->x[d] * cell_N[d] / length[d]) % cell_N[d];

        //      // DEBUG
        // std::cout << tmp[d] << "\t";

    }

    //    //DEBUG FOR-LOOP
    //    std::cout << std::endl;
    //    for (int d=0; d<DIM; d++)
    //        std::cout << "x[" << d << "]: " << i->x[d] << "\t";

    //    std::cout << std::endl << "Corresponding Index: " << J(tmp,cell_N) << std::endl << std::endl;
    return J(tmp,cell_N);

}

real WorldLC::calcBeta(int dimension)
{
    real tmp = 0.0;
    for (std::vector<Cell>::iterator cell = cells.begin(); cell != cells.end (); cell++)
        for (std::list<Particle>::iterator i = cell->particles.begin (); i != cell->particles.end (); i++)
            tmp += sqr(i->v[dimension]);
    return sqrt(thermo_target_temp * (nParticles-1) / (24*tmp));
}

std::ostream& operator << (std::ostream& os, WorldLC& W) 
{
    // Get out some information about the world
    os << W.name << " Dim=" << DIM << " t=" << W.t << " delta_t=" << W.delta_t << " t_end=" << W.t_end
       << " Number of Cells=" << W.cells.size() << " cell_r_cut=" << W.cell_r_cut << std::endl;
    // roll over each Cell
    for (std::vector<Cell>::iterator i = W.cells.begin(); i < W.cells.end(); i++)
    {
        // if there are some particles in this cell, show it off
        if (W.cells[i-W.cells.begin()].particles.size() > 0)
        {
            // we are now in cell# x
            os << "cell[" << i-W.cells.begin() << "] = " ;
            // get out every particle in this cell
            for (std::list<Particle>::iterator j = W.cells[i-W.cells.begin()].particles.begin();j != W.cells[i-W.cells.begin()].particles.end(); j++)
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
        // else observe next cell
        else continue;
    }
    return os;
}

// vim:set et sts=4 ts=4 sw=4 ai ci cin:
