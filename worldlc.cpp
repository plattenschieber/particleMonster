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
    constructParticle(MPI_Particle);
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

    // Calc #cells
    for (int d=0; d<DIM; d++)
    {
        // #cells in dimension = floor(length per cell-cutlength)
        cell_N[d] = (int)(length[d]/cell_r_cut);
        // cell length = world length per #cells
        cell_length[d] = length[d]/cell_N[d];
    }

    // do some parallel stuff

    // Get the number of processes.
    s.numprocs = MPI::COMM_WORLD.Get_size();
    // Get the individual process ID.
    s.myrank = MPI::COMM_WORLD.Get_rank();

    // #procs in dim
    switch(s.numprocs)
    {
        case 1:
            s.N_p[0] = 1;
            s.N_p[1] = 1;
            s.N_p[2] = 1;
            break;
        case 2:
            s.N_p[0] = 2;
            s.N_p[1] = 1;
            s.N_p[2] = 1;
            break;
        case 4:
            s.N_p[0] = 2;
            s.N_p[1] = 2;
            s.N_p[2] = 1;
            break;
        case 6:
            s.N_p[0] = 3;
            s.N_p[1] = 2;
            s.N_p[2] = 1;
            break;
        case 8:
            s.N_p[0] = 2;
            s.N_p[1] = 2;
            s.N_p[2] = 2;
            break;
        default:
            std::cerr << "FAILURE" << std::endl
                      << "You choosed more than 8 processes to work with. Please contact programmer";
            exit(EXIT_FAILURE);
    }

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

    for (int d=0; d<DIM; d++)
    {
        // reset ipTmp
        memcpy(ipTmp, s.ip, sizeof(s.ip));
        // proc isn't at the border
        if (ipTmp[d] > 0)
            ipTmp[d]--;
        // proc is at the border
        else
        {
            ipTmp[d] = NO_NEIGHBOUR;
            if (lower_border[d] == periodic)
                ipTmp[d] = s.N_p[d] - 1;
        }
        // get according number of process
        s.ip_lower[d] = J(ipTmp, s.N_p);

        // do the same stuff for upper border
        memcpy(ipTmp, s.ip, sizeof(s.ip));
        // same procedure here for lower borders
        if (ipTmp[d] < s.N_p[d] - 1)
            ipTmp[d]++;
        else
        {
            ipTmp[d] = NO_NEIGHBOUR;
            if (upper_border[d] == periodic)
                ipTmp[d] = 0;
        }
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
        s.ic_upper_global[d] = s.ic_lower_global[d] - 1  + s.ic_number[d] - 2*s.ic_start[d] ;

    }

    // calc #of all cells incl. bordures
    int nCells = 1;
    for (int d=0; d<DIM; d++)
        nCells *= s.ic_number[d];
    // insert empty cells
    for (int i=0; i<nCells; i++)
        cells.push_back(Cell());

    std::cout << "END OF readParameter()" << std::endl << this;
  //  std::cin >> s.myrank;
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
        if (p.x[d] < 0 || p.x[d] > length[d])
        {
            std::cerr << "<------- FAILURE ------->" << std::endl;
            std::cerr << "Particle left faulty the lower/upper border - in getCellNumber() -> x["
                      << d << "] = " << p.x[d] << std::endl;
            exit(EXIT_FAILURE);
        }
        // compute cell number, AND DON'T FORGET TO ADD IC_START (for the displacement in world)
        tmp[d] = (int) floor(p.x[d] / s.cellh[d]) % s.N_c[d] + s.ic_start[d];
    }
    // return corresponding cell Number
    return J(tmp,s.ic_number);
}


void WorldLC::setCommunication(int dim,
                            int *lower_ic_start, int *lower_ic_stop, int *lower_ic_startreceive, int *lower_ic_stopreceive,
                            int *upper_ic_start, int *upper_ic_stop, int *upper_ic_startreceive, int *upper_ic_stopreceive)
{
    for (int d=0; d<DIM; d++)
    {
        // only bordure
        if (d==dim)
        {
            lower_ic_start[d] = s.ic_start[d];
            lower_ic_stop[d] = 2*lower_ic_start[d];
            lower_ic_startreceive[d] = 0;
            lower_ic_stopreceive[d] = lower_ic_start[d];

            upper_ic_stop[d] = s.ic_stop[d];
            upper_ic_start[d] = s.ic_stop[d] - s.ic_start[d];
            upper_ic_stopreceive[d] = s.ic_start[d] + s.ic_stop[d];
            upper_ic_startreceive[d] = s.ic_stop[d];
        }
        // bordure inclusive
        else if (d>dim)
        {
            lower_ic_startreceive[d] = lower_ic_start[d] = upper_ic_startreceive[d] = upper_ic_start[d] = 0;
            lower_ic_stopreceive[d] = lower_ic_stop[d] = upper_ic_stopreceive[d] = upper_ic_stop[d] = s.ic_start[d] + s.ic_stop[d];
        }
        // w/o bordure
        else
        {
            lower_ic_start[d] = lower_ic_startreceive[d] = upper_ic_start[d] = upper_ic_startreceive[d] = s.ic_start[d];
            lower_ic_stop[d] = lower_ic_stopreceive[d] = upper_ic_stop[d] = upper_ic_stopreceive[d] = s.ic_stop[d];
        }

    }

}



void WorldLC::communicate (bool isForward)
{
    int lower_ic_start[DIM], lower_ic_stop[DIM],
        upper_ic_start[DIM], upper_ic_stop[DIM],
        lower_ic_startreceive[DIM], lower_ic_stopreceive[DIM],
        upper_ic_startreceive[DIM], upper_ic_stopreceive[DIM];

    for (int d=(isForward)?DIM-1:0; (isForward)?d>=0:d<DIM; (isForward)?d--:d++)
    {
        if (isForward)
            setCommunication(d,
                             lower_ic_start, lower_ic_stop, lower_ic_startreceive, lower_ic_stopreceive,
                             upper_ic_start, upper_ic_stop, upper_ic_startreceive, upper_ic_stopreceive);
        else if (!isForward)
            setCommunication(d,
                             lower_ic_startreceive, lower_ic_stopreceive, lower_ic_start, lower_ic_stop,
                             upper_ic_startreceive, upper_ic_stopreceive, upper_ic_start, upper_ic_stop);

        sendReceive (   s.ip_lower[d], lower_ic_start, lower_ic_stop, lower_ic_startreceive, lower_ic_stopreceive,
                        s.ip_upper[d], upper_ic_start, upper_ic_stop, upper_ic_startreceive, upper_ic_stopreceive);

        sendReceive (   s.ip_upper[d], upper_ic_start, upper_ic_stop, upper_ic_startreceive, upper_ic_stopreceive,
                        s.ip_lower[d], lower_ic_start, lower_ic_stop, lower_ic_startreceive, lower_ic_stopreceive);
    }
}

void WorldLC::sendReceive( int lower_proc, int *lower_ic_start,  int *lower_ic_stop, int *lower_ic_startreceive, int* lower_ic_stopreceive,
                           int upper_proc, int *upper_ic_start,  int *upper_ic_stop, int *upper_ic_startreceive, int *upper_ic_stopreceive)
{
    // send and receive infos
    MPI::Status status;
    MPI::Request request;
    // number of particles to be send/received
    int sum_lengthsend = 0, sum_lengthreceive = 0;
    // iterator
    int nCellsSend = 1, k = 0;
    int itCell[DIM];
    std::vector<int> ic_lengthsend, ic_lengthreceive;
    std::vector<Particle> ip_particlesend, ip_particlereceive;

    // both neighbours are there
    if( lower_proc != NO_NEIGHBOUR && upper_proc != NO_NEIGHBOUR )
    {
        // send and receive to/from lowerproc
        for (int d=0; d<DIM; d++)
            nCellsSend *= lower_ic_stop[d] - lower_ic_start[d];
        ic_lengthsend.resize (nCellsSend);
        ic_lengthreceive.resize (nCellsSend);

        k=0;
        int debug;
        Iterate(itCell, lower_ic_start, lower_ic_stop)
        {
            debug = J(itCell,s.ic_number);
            ic_lengthsend[k] = cells[J(itCell,s.ic_number)].particles.size ();
            sum_lengthsend += ic_lengthsend[k++];
        }
        request = MPI::COMM_WORLD.Isend (&(ic_lengthsend.front ()), nCellsSend, MPI::INT, lower_proc, 1);
        MPI::COMM_WORLD.Recv (&(ic_lengthreceive.front ()), nCellsSend, MPI::INT, upper_proc, 1, status);
        request.Wait(status);

        // free the lengthsend for new round
        ic_lengthsend.clear ();
        for (int i=0; i<nCellsSend; i++)
            sum_lengthreceive += ic_lengthreceive[i];

        ip_particlesend.resize (sum_lengthsend);
        ip_particlereceive.resize (sum_lengthreceive);

        k = 0;
        Iterate(itCell, lower_ic_start, lower_ic_stop)
        {
            for (std::list<Particle>::iterator p = cells[J(itCell,s.ic_number)].particles.begin(); p != cells[J(itCell, s.ic_number)].particles.end(); p++)
                ip_particlesend[k++] = *p;
        }

        request = MPI::COMM_WORLD.Isend (&(ip_particlesend.front ()), sum_lengthsend, MPI_Particle, lower_proc, 2);
        MPI::COMM_WORLD.Recv (&(ip_particlereceive.front ()), sum_lengthreceive, MPI_Particle, upper_proc, 2, status);
        request.Wait(status);

        // free the sended particle for new round
        ip_particlesend.clear ();

        // push ic_lengthreceive[kreceive] particle into belonging cell
        k = 0;
        // push all received particles into their according cell
        Iterate( itCell, upper_ic_startreceive, upper_ic_stopreceive )
        {
            for (int i=0; i<ic_lengthreceive[k]; i++)
            {
                Particle *p = new Particle;
                *p = ip_particlereceive[i];
                cells[J(itCell,s.ic_number)].particles.push_back(*p);

                // if belonging cell is inside world (that means, a particle moved from bordure into world)
                int tmp=0;
                for (int d=0; d<DIM; d++)
                {
                    if (itCell[d] > s.ic_start[d] && itCell[d] < s.ic_stop[d])
                        tmp++;
                    else // cell is on bordure, we are in compF case
                        break;
                }
                // don't forget to update particle size
                if (tmp==DIM)
                    nParticles++;
            }
            k++;
        }
        // we received all particles, now clear them
        ic_lengthreceive.clear ();
        ip_particlereceive.clear ();
    }
    // lower neighbour is missing
    else if (upper_proc != NO_NEIGHBOUR )
    {
        // calc number of cells to receive
        for (int d=0; d<DIM; d++)
            nCellsSend *= lower_ic_stop[d] - lower_ic_start[d];
        ic_lengthreceive.resize (nCellsSend);

        // receive displacement
        MPI::COMM_WORLD.Recv (&(ic_lengthreceive.front ()), nCellsSend, MPI::INT, upper_proc, 1);
        for (int i=0; i<nCellsSend; i++)
            sum_lengthreceive += ic_lengthreceive[i];
        ip_particlereceive.resize (sum_lengthreceive);

        // receive particles
        MPI::COMM_WORLD.Recv (&(ip_particlereceive.front ()), sum_lengthreceive, MPI_Particle, upper_proc, 2);
        k = 0;
        Iterate (itCell, upper_ic_startreceive, upper_ic_stopreceive)
        {
            for (int i=0; i<ic_lengthreceive[k]; i++)
            {
                Particle *p = new Particle;
                *p = ip_particlereceive[i];
                cells[J(itCell,s.ic_number)].particles.push_back(*p);
                // if belonging cell is inside world (that means, a particle moved from bordure into world)
                int tmp=0;
                for (int d=0; d<DIM; d++)
                {
                    if (itCell[d] > s.ic_start[d] && itCell[d] < s.ic_stop[d])
                        tmp++;
                    else // cell is on bordure, we are in compF case
                        break;
                }
                // don't forget to update particle size
                if (tmp==DIM)
                    nParticles++; // don't forget to update particle size
            }
            k++;
        }
       ic_lengthreceive.clear ();
       ip_particlereceive.clear ();
    }
    // upper neighbour is missing
    else if (lower_proc != NO_NEIGHBOUR)
    {
        // compute the number of cells in each dim
        for (int d=0; d<DIM; d++)
            nCellsSend *= lower_ic_stop[d] - lower_ic_start[d];
        ic_lengthsend.resize (nCellsSend);

        k = 0;
        // and save number of particles (sum_lengthsend) for each to be sended cell into ic_length
        Iterate (itCell, lower_ic_start, lower_ic_stop)
        {
            ic_lengthsend[k] = cells[J(itCell,s.ic_number)].particles.size();
            sum_lengthsend += ic_lengthsend[k++];
        }

        // send displacement of arriving particles
        MPI::COMM_WORLD.Isend (&(ic_lengthsend.front ()), nCellsSend, MPI::INT, lower_proc, 1);
        ic_lengthsend.clear ();

        // save to be send particles into ip_particlesend
        ip_particlesend.resize (sum_lengthsend);
        k = 0;
        Iterate (itCell, lower_ic_start, lower_ic_stop )
                for (std::list<Particle>::iterator p = cells[J(itCell,s.ic_number)].particles.begin(); p != cells[J(itCell,s.ic_number)].particles.end(); p++)
                ip_particlesend[k++] = *p;
        // and send them to the lower process
        MPI::COMM_WORLD.Isend (&(ip_particlesend.front ()), sum_lengthsend, MPI_Particle, lower_proc, 2);
        ip_particlesend.clear ();
    }
}
void WorldLC::deleteBorderParticles ()
{
    int n[DIM], lower[DIM], upper[DIM];
    for (int i=0; i<DIM; i++)
    {
        lower[i] = 0;
        upper[i] = s.ic_start[i];
        for (int j=0; j<DIM; j++)
            if (i != j)
            {
                lower[j] = 0;
                upper[j] = s.ic_number[j];
            }
        Iterate (n, lower, upper)
            cells[J(n, s.ic_number)].particles.clear();

        lower[i] = s.ic_stop[i];
        upper[i] = s.ic_number[i];
        for (int j=0; j<DIM; j++)
            if (i != j)
            {
                lower[j] = 0;
                upper[j] = s.ic_number[j];
            }
        Iterate (n, lower, upper)
            cells[J(n, s.ic_number)].particles.clear();
    }
}

void WorldLC::constructParticle(MPI::Datatype& MPI_Particle)
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
    // Get out some information about the world and it's SubDomain
    os << W.name << " Dim: " << DIM << " t: " << W.t << " delta_t: " << W.delta_t << " t_end: " << W.t_end
       << " cell_r_cut: " << W.cell_r_cut << std::endl;
    os << "Length: ";
    for (int d=0; d<DIM; d++)
        os << W.length[d] << " ";
    os << std::endl << "Upper borders: ";
    for (int d=0; d<DIM; d++)
        os << W.upper_border[d] << " ";

    os << std::endl << "Lower borders: ";
    for (int d=0; d<DIM; d++)
        os << W.lower_border[d] << " ";
    os << std::endl << "Border types: 0 - leaving, 1 - periodic, 2 - unknown" << std::endl;

    os << std::endl << "Present particles in all " << W.cells.size()
       << " Cells. (#Cells with bordure in " ;
    for (int d=0; d<DIM; d++) os <<  "  dim " << d << ": " << W.s.ic_number[d];
    os << ")" << std::endl;

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
    os << std::endl << "My Rank: " << W.s.myrank << ", position in process Matrix: ";
    for (int d=0; d<DIM; d++) os << W.s.ip[d] << " ";

    os << std::endl << "My lower neighbours (ID's in dth dim): ";
    for (int d=0; d<DIM; d++) os << W.s.ip_lower[d] << " ";

    os << std::endl << "My upper neighbours (ID's in dth dim): ";
    for (int d=0; d<DIM; d++) os << W.s.ip_upper[d] << " ";

    os << std::endl << "My lower global position (of cell): ";
    for (int d=0; d<DIM; d++) os << W.s.ic_lower_global[d] << " ";

    os << std::endl << "My upper global position (of cell): ";
    for (int d=0; d<DIM; d++) os << W.s.ic_upper_global[d] << " ";

    return os;
}

// vim:set et sts=4 ts=4 sw=4 ai ci cin:
