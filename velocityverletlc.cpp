#include "velocityverletlc.hpp"
#include <math.h> 
#include "defines.hpp"

VelocityVerletLC::VelocityVerletLC(WorldLC& _W, LJPotential& _Pot, ObserverXYZ& _O) : VelocityVerlet(_W,_Pot,_O), W(_W), Pot(_Pot) 
{
    // initialize your own World, otherwise implicit cast to World will force us to explicit cast the World every time we use it, to WorldLC
}

VelocityVerletLC::VelocityVerletLC(WorldLC& _W, LJPotential* _Pot, ObserverXYZ& _O) : VelocityVerlet(_W,(*_Pot),_O), W(_W), Pot(*_Pot)
{
    // initialize your own World, otherwise implicit cast to World will force us to explicit cast the World every time we use it, to WorldLC
}

void VelocityVerletLC::compF()
{
    // DEBUG
    real time = W.t;
    // Cell and neighbour cell indices
    int jCell[DIM], nbCell[DIM];
    // check the distance and throw out all that's more far away than rcut
    real dist = 0.0;
    // there is no e_pot in the beginning
    W.e_pot = 0.0;
    // roll over each cell
    for (jCell[0]=0; jCell[0]<W.cell_N[0]; jCell[0]++)
    {
        for (jCell[1]=0; jCell[1]<W.cell_N[1]; jCell[1]++)
        {
            for (jCell[2]=0; jCell[2]<W.cell_N[2]; jCell[2]++)
            {
                // we compute the e_pot for each pair of particles in it's cell including the neighbour cells and add it to the worlds' e_pot...
                // roll over every particle i in actual cell
                for (std::list<Particle>::iterator i = W.cells[J(jCell,W.cell_N)].particles.begin(); i != W.cells[J(jCell,W.cell_N)].particles.end(); i++)
                {
                    // DEBUG:
                    real iTmp[DIM];
                    memcpy(iTmp, i->x, sizeof(i->x));

                    // roll over every neighbour cell
                    for (nbCell[0]=jCell[0]-1; nbCell[0]<=jCell[0]+1; nbCell[0]++)
                    {
                        for (nbCell[1]=jCell[1]-1; nbCell[1]<=jCell[1]+1; nbCell[1]++)
                        {
                            for (nbCell[2]=jCell[2]-1; nbCell[2]<=jCell[2]+1; nbCell[2]++)
                            {
                                //
                                bool leftWorld = false, twistForce = false;
                                bool periodic[DIM] = {false, false, false};

                                // resolve neighbours real position, especially in periodic case
                                int nbTmpCell[DIM];
                                // copy therefor neighbour cells position to an temporay array
                                memcpy(nbTmpCell, nbCell, sizeof(nbCell));

                                // set neighbour to its new place if world is periodic and observed cell is located at the border
                                // or ignore this neighbour cell (leftWorld = true)
                                for (int d=0; d<DIM; d++)
                                {
                                    // accumulate distance between particle and neighbour cell.
                                    // dist += sqr(W.cell_length[d]*nbCell[d] - i->x[d]);

                                    // PERIODIC
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

                                    // LEAVING
                                    else if (nbCell[d]<0 && W.lower_border[d]==W.leaving) leftWorld = true;
                                    else if (nbCell[d]>=W.cell_N[d] && W.upper_border[d]==W.leaving) leftWorld = true;
                                }

                                // compute only if the neighbour cell is inside the world respectivly the world is periodic
                                if(!leftWorld)
                                {
                                    // DEBUG
                                    int watchCell = J(nbTmpCell, W.cell_N);
                                    // foreach particle j in temporary! neighbourcell, compute force
                                    for (std::list<Particle>::iterator j = W.cells[J(nbTmpCell,W.cell_N)].particles.begin(); j != W.cells[J(nbTmpCell,W.cell_N)].particles.end(); j++)
                                    {
                                        // ...except of the computation with itself (i!=j)
                                        if (i!=j)
                                        {
                                            // DEBUG:
                                            // save i's and j's coordinates
                                            real jTmp[DIM];
                                            // copy coordinates
                                            memcpy(jTmp, j->x, sizeof(j->x));
                                            // don't forget to reset the distance
                                            dist = 0.0;
                                            for (int d=0; d<DIM; d++)
                                            {
                                                // IN PERIODIC:
                                                if (periodic[d] && W.cells.size() > 1)
                                                {
                                                    // if nbCell left upper border -> j.x[d] < i.x[d]
                                                    if (j->x[d] < i->x[d])
                                                    {
                                                        // add distance from i to upB and from lowB to j
                                                        dist += sqr((W.length[d]-i->x[d]) + j->x[d]);
                                                    }
                                                    // else nbCell left lower border -> j.x[d] > i.x[d]
                                                    else
                                                    {
                                                        // add distance from lowB to i and from j to upB
                                                        dist += sqr(i->x[d] + (W.length[d]-j->x[d]));
                                                    }

                                                }
                                                // else cell is only cell in periodic case or we are in the world and the distance is calculated as usual
                                                // TODO: case Periodic and 1 cell!! -> particle1 on the right, particle2 on the left -> dist periodic is shorter than normal dist
                                                else
                                                    dist += sqr(j->x[d] - i->x[d]);
                                            }
                                            // DEBUG
                                            std::cout << "distance between i and j: " << dist << std::endl;
                                            // only particles which are closer than rcut
                                            if (dist <= W.cell_r_cut)
                                                // computes the force between particle i and j and add it to our potential
                                                W.e_pot += Pot.force(*i, *j, dist, W.epsilon, W.sigma);
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    
}

void VelocityVerletLC::updateV()
{
    // there is no e_kin in the beginning
    W.e_kin = 0.0;
    // roll over every cell
    for (std::vector<Cell>::iterator cell =  W.cells.begin(); cell < W.cells.end(); cell++)
    {
        // foreach cell go through it's particles...
        for (std::list<Particle>::iterator i = cell->particles.begin(); i != cell->particles.end(); i++)
        {
            // ...and over every dimension of particle i
            for (unsigned int d=0; d<DIM; d++)
            {
                // compute new velocity in dimension d
                i->v[d] += .5*(i->F_old[d] + i->F[d])*W.delta_t/i->m;
                // if we want to check the temperatur regulary
                if (fmod(W.t,W.thermo_step_interval) == 0 && W.isThermoStartTemp)
                    // multiply velocity by beta
                    i->v[d] *= W.calcBeta();
                // add now the pro rata e_kin
                W.e_kin += .5*i->m*sqr(i->v[d]);
            }
        }
    }


}

void VelocityVerletLC::updateX()
{
    // if the flag is checked, push the particle in the last round into it's new position
    bool doIt = false;
    bool innerWorld = true;

    // roll over every cell
    real time = W.t;
    for (std::vector<Cell>::iterator cell =  W.cells.begin(); cell < W.cells.end(); cell++)
    {
        // foreach cell go through it's particles...
        for (std::list<Particle>::iterator i = cell->particles.begin(); i != cell->particles.end(); i++)
        {
            Particle &p = *i;
            // DEBUG at first get out every particle and it's cell number
            std::cout << " Cell[" << W.getCellNumber(i) << "]"
                      << ".particle["  <<  i->ID  << "]";
            std::cout << "x -> ";
            for (int d=0; d<DIM; d++) std::cout << i->x[d] << " ";
            std::cout << "v -> ";
            for (int d=0; d<DIM; d++) std::cout << i->v[d] << " ";
            std::cout << "F -> ";
            for (int d=0; d<DIM; d++) std::cout << i->F[d] << " ";
            std::cout << std::endl;


            // if the flag is checked, push the particle in the last round into it's new position
            doIt = false;
            // if flag is not checked, particle is we are at the border
            innerWorld = true;


            // compare Cell Numbers of (maybe) moved particle i
            int checkCell = W.getCellNumber(i);
            // 	..at first calc new position in every dimension
            for (unsigned int d=0; d<DIM; d++)
            {
                // computing new location of the particle i if it's leaving the world, elsewise just call handleBorders (-lc version) in the end

                i->x[d] += W.delta_t*i->v[d] + (.5*i->F[d]*sqr(W.delta_t)) / i->m;


                //                std::cout << "Cell[" << W.getCellNumber(i) << "]"
                //                          << ".particle["  <<  i->ID  << "]"
                //                          << ".x[" << d << "]=" << i->x[d] << std::endl;

                // save last force...
                i->F_old[d] = i->F[d];
                // ... and don't forget to set the actual force to zero
                i->F[d] = 0;
            }
            //std::cout << W.t << " Cell[" << W.getCellNumber(i) << "]"
            //<< ".particle["  <<  i->ID  << "]" << std::endl;
            // then check if particle left its cell and handle moving issues (respective border issues)
            if (W.getCellNumber(i) != checkCell)
            {
                // check it in every dimension
                for (unsigned int d = 0; d<DIM; d++)
                {
                    // periodic - position = position % worldlength
                    if (i->x[d] >= W.length[d] && W.upper_border[d] == W.periodic)
                    {
                        // DEBUG:
                        std::cout << "New position (oben Raus Untenwiederrein)" << std::endl;
                        // new position is at the beginning of world plus the overhead that x left the world
                        i->x[d] = fmod(i->x[d], W.length[d]);

                        // there is a particle which left the world at one of its sides
                        doIt = true;
                    }
                    // periodic - position = position % worldlength
                    else if (i->x[d] < 0 && W.lower_border[d] == W.periodic)
                    {
                        // DEBUG:
                        std::cout << "New position (unten Raus Obenwiederrein)" << std::endl;
                        // new position is end of world minus overhead that x left the world
                        i->x[d] =  W.length[d] - fabs(fmod(i->x[d], W.length[d]));

                        // there is a particle which left the world at one of its sides
                        doIt = true;
                    }
                    // leaving - it just bumps out
                    else if (i->x[d] >= W.length[d] && W.upper_border[d] == W.leaving)
                    {
                        // DEBUG:
                        std::cout << "New position (oben raus Wegvomfenster): " << std::endl;
                        // regardless in which dimension, just erase
                        i = cell->particles.erase(i);
                        i--;
                        W.nParticles--;
                        // don't forget to set dimension to DIM  or you will handle another particle in the wrong dimension
                        break;
                    }
                    // leaving - it just bumps out
                    else if (i->x[d] < 0  && W.lower_border[d] == W.leaving)
                    {
                        // DEBUG:
                        std::cout << "New position (unten raus Wegvomfenster): " << std::endl;
                        // regardless in which dimension, just erase
                        i = cell->particles.erase(i);
                        i--;
                        W.nParticles--;
                        break;
                    }

                    // calc the position in all dimension and THEN push it into the right cell
                    if (d==DIM-1 && doIt==true)
                    {
                        // DEBUG:
                        //std::cout << "New position (oben raus Periodisch): " << i->x[d] << std::endl;

                        // Add particle i to its corresponding cell
                        W.cells[W.getCellNumber(i)].particles.push_back(*i);
                        i = cell->particles.erase(i);
                        i--;
                        break;
                    }
                    if (d==DIM-1 && doIt==false)
                    {
                        W.particles.push_back(*i);
                        i = cell->particles.erase(i);
                        i--;
                        break;
                    }
                }
            }

        }
    }
    // and now add the particles again to their belonging cells
    for (std::list<Particle>::iterator i = W.particles.begin(); i != W.particles.end(); i++)
    {
        //W.cells[W.getCellNumber(i)].particles.push_back(cell->particles[i-cell->particles.begin()]);
        W.cells[W.getCellNumber(i)].particles.push_back(*i);
        i = W.particles.erase(i);
        i--;
    }

}
// vim:set et sts=4 ts=4 sw=4 ai ci cin:
