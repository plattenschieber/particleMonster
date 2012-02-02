#include "velocityverletlc.hpp"
#include <math.h> 
#include "defines.hpp"
#include <sstream>


VelocityVerletLC::VelocityVerletLC(WorldLC& _W, LJPotential& _Pot, ObserverXYZ& _O) : VelocityVerlet(_W,_Pot,_O), W(_W), Pot(_Pot)
{
    // initialize your own World, otherwise implicit cast to World will force us to explicit cast the World every time we use it, to WorldLC
}
VelocityVerletLC::VelocityVerletLC(WorldLC& _W, LJPotential* _Pot, ObserverXYZ& _O) : VelocityVerlet(_W,(*_Pot),_O), W(_W), Pot(*_Pot)
{
    // initialize your own World, otherwise implicit cast to World will force us to explicit cast the World every time we use it, to WorldLC
}

void VelocityVerletLC::simulate ()
{
    // initialize forces
    //initF();
    // watch out for data in before first step
    //O.notify ();
    // call base class method
    VelocityVerlet::simulate ();
    // call the destructor explicitly
    O.~ObserverXYZ ();
}

void VelocityVerletLC::initF()
{
    for (std::vector<Cell>::iterator cell =  W.cells.begin(); cell < W.cells.end(); cell++)
        if (!cell->particles.empty())
            // foreach cell go through it's particles...
            for (std::list<Particle>::iterator i = cell->particles.begin(); i != cell->particles.end(); i++)
                for (int d=0; d<DIM; d++)
                {
                    i->F_old[d] = i->F[d];
                    i->F[d] = 0;
                }
    W.e_pot = 0;
}

void VelocityVerletLC::compF()
{
    // fill bordure with particles
    W.communicate (true);

    // Cell and neighbour cell indices
    int jCell[DIM], nbCell[DIM];

    // there is no e_pot in the beginning
    W.e_pot = 0.0;

    // roll over each cell
    Iterate (jCell, W.s.ic_start, W.s.ic_stop)
    {
        // we compute the e_pot for each pair of particles in it's cell including the neighbour cells and add it to the worlds' e_pot...
        // roll over every particle i in actual cell
        for (std::list<Particle>::iterator i = W.cells[J(jCell,W.s.ic_number)].particles.begin(); i != W.cells[J(jCell,W.s.ic_number)].particles.end(); i++)
        {
            // roll over every neighbour cell
            Iterate (nbCell, jCell -1, jCell +1)
            {
                // mark that a neighbour is outside the world
                bool leftWorld = false;
                // mark that a neighbour left the  world at a specific periodic border
                bool periodic[DIM];
                for (int d=0; d<DIM; d++)
                    periodic[d] = false;

                // resolve neighbours real position, especially in periodic case
                int nbTmpCell[DIM];
                //copy therefor neighbour cells position to an temporay array
                memcpy(nbTmpCell, nbCell, sizeof(nbCell));

                // set neighbour to its new place if world is periodic and observed cell is located at the border
                // or ignore this neighbour cell (leftWorld = true)
                for (int d=0; d<DIM; d++)
                {
                    // PERIODIC CASES:
                    if (nbCell[d]<0 && W.lower_border[d]==W.periodic)
                    {
                        nbTmpCell[d] = W.cell_N[d]-1;
                        periodic[d] = true;
                    }
                    else if (nbCell[d]>=W.cell_N[d] && W.upper_border[d]==W.periodic)
                    {
                        nbTmpCell[d]=0;
                        periodic[d] = true;
                    }

                    // LEAVING CASES:
                    else if (nbCell[d]<0 && W.lower_border[d]==W.leaving) leftWorld = true;
                    else if (nbCell[d]>=W.cell_N[d] && W.upper_border[d]==W.leaving) leftWorld = true;
                    else if (nbCell[d]<0 && W.lower_border[d]==W.unknown)
                    {
                        std::cerr << "<------- FAILURE ------->" << std::endl;
                        std::cerr << "Please specify lower border in Dimension: " << d+1;
                        exit(EXIT_FAILURE);
                    }
                    else if (nbCell[d]>=W.cell_N[d] && W.upper_border[d]==W.unknown)
                    {
                        std::cerr << "<------- FAILURE ------->" << std::endl;
                        std::cerr << "Please specify upper border in Dimension: " << d+1;
                        exit(EXIT_FAILURE);
                    }
                }

                // compute only if the neighbour cell is inside the world
                if(!leftWorld)
                {
                    // foreach particle j in temporary! neighbourcell compute force
                    for (std::list<Particle>::iterator j = W.cells[J(nbTmpCell,W.s.ic_number)].particles.begin(); j != W.cells[J(nbTmpCell,W.s.ic_number)].particles.end(); j++)
                    {
                        // ...except of the computation with itself (i!=j)
                        if (i!=j)
                        {
                            // check distance and drop all particle which are farther than rcut
                            real dist = 0.0;
                            // save a direction vector for the distance
                            real dirV[DIM];
                            for (int d=0; d<DIM; d++)
                                dirV[d] = 0.0;

                            // compute distance and directional vector of i and j
                            for (int d=0; d<DIM; d++)
                            {
                                // IN PERIODIC:
                                if (periodic[d] && W.cells.size() > 1)
                                {
                                    // if nbCell left upper border -> j.x[d] < i.x[d]
                                    if (j->x[d] < i->x[d])
                                        // add distance from i to upB and from lowB to j
                                        dirV[d] = (W.length[d] - i->x[d]) + j->x[d];

                                    // else nbCell left lower border -> j.x[d] > i.x[d]
                                    else
                                        // add distance from lowB to i and from j to upB (times -1, because of the over border handling)
                                        dirV[d] = -1*(i->x[d] + (W.length[d] - j->x[d]));
                                }
                                // PERIODIC 1 cell:
                                else if (periodic[d] && W.cells.size () == 1)
                                {
                                    if( (j->x[d] - i->x[d]) > 0.5*W.s.cellh[d])
                                        // and update direction vector
                                        dirV[d] = i->x[d] - j->x[d];
                                    else
                                        dirV[d] = j->x[d] - i->x[d];
                                }
                                // IN NONPERIODIC:
                                else
                                {
                                    // and update direction vector
                                    dirV[d] = j->x[d] - i->x[d];
                                }

                                // And now, accumulate the distance...
                                dist += sqr(dirV[d]);
                            }

                            // only particles which are closer than rcut, flow into the computation
                            if (dist <= sqr(W.cell_r_cut))
                                // computes the force between particle i and j and add it to our potential
                                W.e_pot += Pot.force(*i, *j, dist, dirV, W.epsilon, W.sigma);
                        }
                    }
                } // end if (inside world)
            } //end neighbour cell Iterate
        } //end iterator over all particles i
    } //end Iterate(jCell)

    // and delete the bordure again
    W.deleteBorderParticles ();
    
}

void VelocityVerletLC::updateV()
{

    // there is no e_kin in the beginning
    W.e_kin = 0.0;
    // velocity scaling
    real beta = 1.0;
    // roll over every cell
    for (std::vector<Cell>::iterator cell =  W.cells.begin(); cell < W.cells.end(); cell++)
        // foreach cell go through it's particles...
        for (std::list<Particle>::iterator i = cell->particles.begin(); i != cell->particles.end(); i++)
            // ...and over every dimension of particle i
            for (unsigned int d=0; d<DIM; d++)
            {
                // compute new velocity in dimension d
                i->v[d] += .5*(i->F_old[d] + i->F[d])*W.delta_t/i->m;
                // add now the pro rata e_kin
                W.e_kin += .5*i->m*sqr(i->v[d]);

                // VELOCITY SCALING
                if (W.isThermoStartTemp && (W.step % W.T_Step == 0) && (fabs(W.T_D - W.T) > 10e-6) )
                {
                    beta = sqrt(W.T_D * (W.nParticles-1) / (48*W.e_kin));
                    std::cout << "VELOTCITY SCALING" << std::endl;
                    // multiply velocity by beta
                    i->v[d] *= beta;
                    // and scale kinetc energy
                    W.e_kin *= beta;
                }
            } // END d<DIM loop
    // compute total energy
    W.e_tot = W.e_kin + W.e_pot;
    // set new temperature
    W.T *= sqr(beta);
}

void VelocityVerletLC::updateX()
{
    // if the flag is checked, push the particle in the last round into it's new position
    bool doIt = false;
    bool innerWorld = true;
    int jCell[DIM], tmpCell[DIM];
    std::vector<int> kInsert;
    // roll over every cell
    Iterate(jCell, W.s.ic_start, W.s.ic_stop)
    {
        // foreach cell go through it's particles...
        for (std::list<Particle>::iterator i = W.cells[J(jCell, W.s.ic_number)].particles.begin(); i != W.cells[J(jCell, W.s.ic_number)].particles.end(); i++)
        {
            // if the flag is checked, push the particle in the last round into it's new position
            doIt = false;
            // if flag is not checked, particle is we are at the border
            innerWorld = true;
            // compare Cell Numbers of (maybe) moved particle i
            int checkCell = W.getCellNumber(*i);
            // 	..at first calc new position in every dimension
            for (unsigned int d=0; d<DIM; d++)
            {
                // computing new location of the particle i if it's leaving the world, elsewise just call handleBorders (-lc version) in the end
                i->x[d] += W.delta_t*i->v[d] + (.5*i->F[d]*sqr(W.delta_t)) / i->m;
                // save last force...
                i->F_old[d] = i->F[d];
                // ... and don't forget to set the actual force to zero
                i->F[d] = 0;
            }

            // if particle changed cell
            if (W.getCellNumber(*i) != checkCell)
            {
                // check it in every dimension
                for (unsigned int d = 0; d<DIM; d++)
                {
                    // leaving - it just bumps out
                    if (i->x[d] > W.length[d] && W.upper_border[d] == W.leaving)
                    {
                        // regardless in which dimension, just erase
                        i = W.cells[J(jCell, W.s.ic_number)].particles.erase(i);
                        i--;
                        W.nParticles--;
                        // don't forget to break, or you will handle another particle in the wrong dimension
                        break;
                    }
                    // leaving - it just bumps out
                    else if (i->x[d] < 0  && W.lower_border[d] == W.leaving)
                    {
                        // regardless in which dimension, just erase
                        i = W.cells[J(jCell, W.s.ic_number)].particles.erase(i);
                        i--;
                        W.nParticles--;
                        break;
                    }
                    // unknown - it just bumps out
                    else if (i->x[d] > W.length[d] && W.upper_border[d] == W.unknown)
                    {
                        std::cerr << "<------- FAILURE ------->" << std::endl;
                        std::cerr << "Please specify upper border in Dimension: " << d;
                        exit(EXIT_FAILURE);
                    }
                    // unknown - it just bumps out
                    else if (i->x[d] < 0  && W.lower_border[d] == W.unknown)
                    {
                        std::cerr << "<------- FAILURE ------->" << std::endl;
                        std::cerr << "Please specify lower border in Dimension: " << d;
                        exit(EXIT_FAILURE);
                    }
                }
                // push particles into waiting pot
                int debug = W.getCellNumber (*i);
                debug = J(jCell, W.s.ic_number);
                W.particles.push_back(*i);
                i = W.cells[J(jCell, W.s.ic_number)].particles.erase(i);
                i--;
                break;
            }

        }
    }

    // calc according cell
    int inCell[DIM];
    // is particle still in our Subdomain?
    bool isInSubdomain = true;
    // add particles to new cells. If went onto bordure, decrease nParticles
    for (std::list<Particle>::iterator i = W.particles.begin(); i != W.particles.end(); i++)
    {

        // check if particle doesn't belong to our subdomain
        for (int d=0; d<DIM; d++)
        {
            // calc local cell number and compare to local indices of subdomain
            inCell[d] = (int)floor( i->x[d]/W.s.cellh[d] ) - W.s.ic_lower_global[d] + W.s.ic_start[d];
            if( inCell[d] < W.s.ic_start[d] || inCell[d] > W.s.ic_stop[d] )
                isInSubdomain = false;
        }
        // when not in inner SubDomain, there was a particle gone
        if (!isInSubdomain)
           W.nParticles--;
        W.cells[J(inCell, W.s.ic_number)].particles.push_back(*i);
        i = W.particles.erase(i);
        i--;
    }

    // backward direction - send bordure to procs
    W.communicate (false);
    // and delete them
    W.deleteBorderParticles ();

}


// vim:set et sts=4 ts=4 sw=4 ai ci cin:
