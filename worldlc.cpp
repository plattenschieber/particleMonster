#include "worldlc.hpp"
#include <mpi.h>
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
    // construct it, so we can use it!
    construct_particle(MPI_Particle);
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

    // do some parallel stuff

    // Get the number of processes.
    s.numprocs = MPI::COMM_WORLD.Get_size();
    // Get the individual process ID.
    s.myrank = MPI::COMM_WORLD.Get_rank();

    // #procs in dim
    s.N_p[0] = s.N_p[1] = 1; s.N_p[2] = 1;

    // get position of actual process in the grid
    Jinv(s.myrank, s.N_p, s.ip);

    // save the global cell displacement of every process
    int *displ[DIM];

    // compute cell number of inner domain in dth dimension
    for (int d=0; d<DIM; d++)
    {
        // we need #procs place for displacements
        displ[d] = (int*)malloc (s.N_p[d] * sizeof(int));
        // s.N_c which needs to be divided to procs
        s.N_c[d] = cell_N[d];

        // calc until you reach your process position
        for(int i=0; i<s.ip[d];i++)
        {
            // save the left over cells as displacement for ith process
            displ[d][i] = cell_N[d] - s.N_c[d];
            // and compute the left over for the next step
            s.N_c[d] -= round(s.N_c[d]/(s.N_p[d]-i));
        }

        // save last left over and continue further down
        int tmp = s.N_c[d];

        // the last calculation equals the number of cells for our process
        s.N_c[d] = round(s.N_c[d]/(s.N_p[d]-s.ip[d]));

        // now compute remaining displacements
        for (int i=s.ip[d]; i<s.N_p[d]; i++)
        {
            displ[d][i] = cell_N[d] - tmp;
            tmp -= round(tmp/(s.N_p[d]-i));
        }
    }

    // temporary placeholder for resolving coordinates of neighbours
    int ipTmp[DIM];
    // save current coordinates in grid to ipTmp
    memcpy(ipTmp, s.ip, sizeof(s.ip));

    for (int d=0; d<DIM; d++)
    {
        // lower neighbour in dimension d (plus s.N_p[d], cause of modulo disability to calc negatives)
        ipTmp[d] = (s.ip[d] - 1 + s.N_p[d]) % s.N_p[d];
        // get according number of process
        s.ip_lower[d] = J(ipTmp, s.N_p);
        // upper neigbour in dimension d
        ipTmp[d] = (s.ip[d] + 1) % s.N_p[d];
        // get according number of process
        s.ip_upper[d] = J(ipTmp, s.N_p);

        // the cells edge length is worlds edge length per #cells
        s.cellh[d] = length[d] / cell_N[d];
        // bordure width - equals the first cell inside subdomain
        s.ic_start[d] = (int) ceil(cell_r_cut / s.cellh[d]);
        // first cell in upper bordure, so you can use it as for i=s.ic_start; i<s.ic_stop
        s.ic_stop[d] = s.ic_start[d] + s.N_c[d];

        // number of cells in dth dim incl. bordure
        // (s.ic_stop[d] - s.ic_start[d]) + 2*s.ic_start[d] is the same as:
        s.ic_number[d] = s.ic_stop[d] + s.ic_start[d];
        // the lower global Index of the first cell, is the dth displacement of the corresponding process
        s.ic_lower_global[d] = displ[d][s.ip[d]];
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
        // getCellNumber(i) gives belonging cellnumber, push into this cell our actual particle i
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
    // temporary array
    int tmp[3] = {0,0,0};
    //    // DEBUG Table
    //    std::cout << "Cell coordinate: " ;
    for (int d=0; d<DIM; d++)
    {
        // handle particle outside the world failure
        if (i->x[d] < 0)
            return -1;
        tmp[d] = (int) floor(i->x[d] / cell_length[d]) % cell_N[d];
    }
    // return corresponding cell Number
    return J(tmp,cell_N);
}

real WorldLC::calcBeta()
{
    return sqrt(thermo_target_temp * (nParticles-1) / (48*e_kin));
}

void WorldLC::SetCommunication(SubDomain *s, int dim,
                            int *lower_ic_start, int *lower_ic_stop, int *lower_ic_startreceive, int *lower_ic_stopreceive,
                            int *upper_ic_start, int *upper_ic_stop, int *upper_ic_startreceive, int *upper_ic_stopreceive)
{
    for (int d=0; d<DIM; d++)
    {
        // only bordure
        if (d==dim)
        {
            lower_ic_start[d] = s->ic_start[d];
            lower_ic_stop[d] = 2*lower_ic_start[d];
            lower_ic_startreceive[d] = 0;
            lower_ic_stopreceive[d] = lower_ic_start[d];

            upper_ic_stop[d] = s->ic_stop[d];
            upper_ic_start[d] = s->ic_stop[d] - s->ic_start[d];
            upper_ic_stopreceive[d] = s->ic_start[d] + s->ic_stop[d];
            upper_ic_startreceive[d] = s->ic_stop[d];
        }
        // bordure inclusive
        else if (d>dim)
        {
            lower_ic_startreceive[d] = lower_ic_start[d] = upper_ic_startreceive[d] = upper_ic_start[d] = 0;
            lower_ic_stopreceive[d] = lower_ic_stop[d] = upper_ic_stopreceive[d] = upper_ic_stop[d] = s->ic_start[d] + s->ic_stop[d];
        }
        // w/o bordure
        else
        {
            lower_ic_start[d] = lower_ic_startreceive[d] = upper_ic_start[d] = upper_ic_startreceive[d] = s->ic_start[d];
            lower_ic_stop[d] = lower_ic_stopreceive[d] = upper_ic_stop[d] = upper_ic_stopreceive[d] = s->ic_stop[d];
        }

    }

}

void WorldLC::Communication (Cell *grid, SubDomain *s, bool isForward)
{
    int lower_ic_start[DIM], lower_ic_stop[DIM],
        upper_ic_start[DIM], upper_ic_stop[DIM],
        lower_ic_startreceive[DIM], lower_ic_stopreceive[DIM],
        upper_ic_startreceive[DIM], upper_ic_stopreceive[DIM];

    for (int d=(isForward)?DIM-1:0; (isForward)?d<DIM:d>=0; (isForward)?d--:d++)
    {
        if(isForward)
            SetCommunication(s, d, lower_ic_start, lower_ic_stop, lower_ic_startreceive, lower_ic_stopreceive,
                              upper_ic_start, upper_ic_stop, upper_ic_startreceive, upper_ic_stopreceive);
        else if(!isForward)
            SetCommunication(s, d, lower_ic_startreceive, lower_ic_stopreceive, lower_ic_start, lower_ic_stop,
                            upper_ic_startreceive, upper_ic_stopreceive, upper_ic_start, upper_ic_stop);

        sendReceive (grid, s->ic_number, s->ip_lower[d], lower_ic_start, lower_ic_stop, lower_ic_startreceive, lower_ic_stopreceive,
                                        s->ip_upper[d], upper_ic_start, upper_ic_stop, upper_ic_startreceive, upper_ic_stopreceive);

    }
}

void WorldLC::sendReceive(Cell *grid, int *ic_number,
                              int lower_proc, int *lower_ic_start,  int *lower_ic_stop, int *lower_ic_startreceive, int* lower_ic_stopreceive,
                              int upper_proc, int *upper_ic_start,  int *upper_ic_stop, int *upper_ic_startreceive, int *upper_ic_stopreceive)
{
    MPI::Status status;
    int sum_lengthsend = 0, sum_lengthreceive = 0;
    int k = 0, kreceive = 0, ncs = 1;
    int *ic_lengthsend = NULL, *ic_lengthreceive = NULL, ic[DIM];
    Particle *ip_particlesend = NULL, *ip_particlereceive = NULL;

    // send and receive to/from lowerproc
    for (int d=0; d<DIM; d++)
        ncs *= lower_ic_stop[d] - lower_ic_start[d];
    ic_lengthsend = (int*)malloc(ncs*sizeof(*ic_lengthreceive));

    //iterate over
    for (int i=0; i<DIM; i++)
    {
        ic_lengthsend[k] = grid[J(ic,ic_number)].particles.size ();
        sum_lengthsend += ic_lengthsend[k++];
    }
    MPI::COMM_WORLD.Isend (ic_lengthsend, ncs, MPI_INT, lower_proc, 1);
    MPI::COMM_WORLD.Recv (ic_lengthreceive, ncs, MPI::INT, upper_proc, 1, status);
    //MPI::Request::Wait();
    //status
    free(ic_lengthsend);
    for (k=0; k<ncs; k++)
        sum_lengthreceive += ic_lengthreceive[k];
    sum_lengthsend *= sizeof(*ip_particlesend);
    ip_particlesend = (Particle*)malloc(sum_lengthsend);
    sum_lengthreceive *= sizeof(*ip_particlereceive);
    ip_particlereceive = (Particle*)malloc (sum_lengthreceive);
    k=0;
    //iterate()
    //for (int i=0; i<DIM; i++)
        //for (std::vector<Cell>::iterator i = grid[J(ic,ic_number)]; i!= grid[J(ic,ic_number)]particles.end (); i++)
            //ip_particlesend[k++] = *i;
    MPI::COMM_WORLD.Isend (ip_particlesend, sum_lengthsend, MPI::CHAR, lower_proc, 2);
    MPI::COMM_WORLD.Recv (ip_particlereceive, sum_lengthreceive, MPI::CHAR, upper_proc, 2, status);
    //MPI::Request::Wait();
    free (ip_particlesend);
    kreceive = k = 0;

    for (int d=0; d<DIM; d++)
    {
        for (int icp=0; icp<ic_lengthreceive[kreceive]; icp++)
        {
            std::vector<Particle> i;
            i.push_back (ip_particlereceive[k++]);
        }
        kreceive++;
    }

    free (ic_lengthreceive);
    free (ip_particlereceive);

}

void WorldLC::construct_particle(MPI::Datatype& MPI_Particle)
{
// Initialize Particle
  Particle p;
  p.ID = 0;
  p.m = 0.0;
  p.x[0] = p.x[1] = p.x[2] = 0.0;
  p.v[0] = p.v[1] = p.v[2] = 0.0;
  p.F[0] = p.F[1] = p.F[2] = 0.0;
  p.F_old[0] = p.F_old[1] = p.F_old[2] = 0.0;

  // build new mpi datatype
  MPI::Datatype type[6] = {MPI::INT, MPI::DOUBLE, MPI::DOUBLE, MPI::DOUBLE, MPI::DOUBLE, MPI::DOUBLE};
  int blocklen[6] = {1,1,DIM,DIM,DIM,DIM};
  MPI::Aint disp[6]; // displacements
  MPI::Aint base;

  /* compute displacements */
  disp[0] = MPI::Get_address(&p);
  disp[1] = MPI::Get_address(&p.m);
  disp[2] = MPI::Get_address(&p.x);
  disp[3] = MPI::Get_address(&p.v);
  disp[4] = MPI::Get_address(&p.F);
  disp[5] = MPI::Get_address(&p.F_old);

  base = disp[0];
  for (unsigned i=0; i<6; i++) disp[i] -= base;

  MPI_Particle = MPI::Datatype::Create_struct(6, blocklen, disp, type);
  MPI_Particle.Commit();
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
    }
    return os;
}

// vim:set et sts=4 ts=4 sw=4 ai ci cin:
