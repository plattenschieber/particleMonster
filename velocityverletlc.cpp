#include "velocityverletlc.hpp"
#include <math.h> 
#include "defines.hpp"



// another macro to expand - this time a for loop based on dimension DIM
#

VelocityVerletLC::VelocityVerletLC(WorldLC& _W, LJPotential& _Pot, ObserverXYZ& _O) : VelocityVerlet(_W,_Pot,_O), W(_W), Pot(_Pot) 
{
    // initialize your own World, otherwise implicit cast to World will force us to explicit cast the World every time we use it, to WorldLC
}
VelocityVerletLC::VelocityVerletLC(WorldLC& _W, LJPotential* _Pot, ObserverXYZ& _O) : VelocityVerlet(_W,(*_Pot),_O), W(_W), Pot(*_Pot)
{
    // initialize your own World, otherwise implicit cast to World will force us to explicit cast the World every time we use it, to WorldLC
}

void VelocityVerletLC::comp_F()
{
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
                for (std::vector<Particle>::iterator i = W.cells[J(jCell,W.cell_N)].particles.begin(); i < W.cells[J(jCell,W.cell_N)].particles.end(); i++)
                {
                    // roll over every neighbour cell
                    for (nbCell[0]=jCell[0]-1; nbCell[0]<=jCell[0]; nbCell[0]++)
                    {
                        for (nbCell[1]=jCell[1]-1; nbCell[1]<=jCell[1]; nbCell[1]++)
                        {
                           for (nbCell[2]=jCell[2]-1; nbCell[2]<=jCell[2]; nbCell[2]++)
                           {
                              bool leftWorld = false;
                              bool periodic[DIM] = {false, false, false};
                              // don't forget to reset the distance
                              //dist = 0.0;
                              // set neighbour to its new place if world is periodic and observed cell is located at the border  
                              for (int d=0; d<DIM; d++)
                              {
                                 // accumulate distance between particle and neighbour cell. 
                                 // dist += sqr(W.cell_length[d]*nbCell[d] - i->x[d]);

                                 // PERIODIC
                                 if (nbCell[d]<0 && W.lower_border[d]==W.periodic)
                                 {
                                     nbCell[d] = W.cell_N[d]-1;
                                     periodic[d] = true;
                                 }
                                 else if (nbCell[d]>=W.cell_N[d] && W.upper_border[d]==W.periodic)
                                 {
                                     nbCell[d]=0;
                                     periodic[d] = true;
                                 }
                                 // LEAVING
                                 if (nbCell[d]<0 && W.lower_border[d]==W.leaving) leftWorld = true;
                                 else if (nbCell[d]>=W.cell_N[d] && W.upper_border[d]==W.leaving) leftWorld = true;


                              }
                              // If neighbour cell is too far away, no calc of inner (more far away) particles are needed! 
                              // if (dist <= W.cell_r_cut)

                              // if the cell isn't outside the world compute!
                              if(!leftWorld)
                              {
                                   // foreach particle j in neighbourcell compute force
                                   for (std::vector<Particle>::iterator j = W.cells[J(nbCell,W.cell_N)].particles.begin(); j < W.cells[J(nbCell,W.cell_N)].particles.end(); j++)
                                   {
                                   // ...except of the computation with itself (i!=j)
                                     if (i!=j)
                                     {
                                         // don't forget to reset the distance
                                         dist = 0.0;
                                         for (int d=0; d<DIM; d++)
                                         {
                                             // IN PERIODIC:
                                             if (periodic[d] && W.cells.size() > 1)
                                             {
                                                // if cell left upper border -> j.x[d] < i.x[d]
                                                if (j->x[d] < i->x[d])
                                                {
                                                     // calc cell coordinate in which i lies
                                                     int tmp = i->x[d] * W.cell_N[d] / W.length[d];
                                                     // add the distance from i to the right border
                                                     dist += (W.cell_length[d] - (i->x[d] - tmp));

                                                     // calc cell coordinate in which j lies
                                                     tmp = j->x[d] * W.cell_N[d] / W.length[d];
                                                     // add the distance from j to the left border
                                                     dist += j->x[d] - tmp;
                                                }
                                                // else if cell left lower border -> j.x[d] > i.x[d]
                                                else /*if (j->x[d] > i->x[d])*/
                                                {
                                                    // calc cell coordinate in which i lies
                                                    int tmp = j->x[d] * W.cell_N[d] / W.length[d];
                                                    // add the distance from i to the right border
                                                    dist += sqr(W.cell_length[d] - (i->x[d] - tmp));

                                                    // calc cell coordinate in which j lies
                                                    tmp = i->x[d] * W.cell_N[d] / W.length[d];
                                                    // add the distance from j to the left border
                                                    dist += sqr(j->x[d] - tmp);
                                                }

                                             }
                                             // else cell is only cell in periodic case or we are in the world and the distance is calculated as usual
                                             else
                                                dist += sqr(j->x[d] - i->x[d]);
                                         }
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

void VelocityVerletLC::update_V()
{
    // there is no e_kin in the beginning
	W.e_kin = 0.0;
	// roll over every cell	
    for (std::vector<Cell>::iterator cell =  W.cells.begin(); cell < W.cells.end(); cell++)
    {
        // foreach cell go through it's particles... 
        for (std::vector<Particle>::iterator i = cell->particles.begin(); i < cell->particles.end(); i++)
        {
            // ...and over every dimension of particle i
            for (unsigned int d=0; d<DIM; d++)
            {
                // compute new velocity in dimension d
                i->v[d] += .5*(i->F_old[d] + i->F[d])*W.delta_t/i->m;
                // add now the pro rata e_kin 
                W.e_kin += .5*i->m*sqr(i->v[d]);
            }
        }
    }


}

void VelocityVerletLC::update_X()
{
    // if the flag is checked, push the particle in the last round into it's new position
    bool doIt = false;
    bool innerWorld = true;
    // roll over every cell
   	for (std::vector<Cell>::iterator cell =  W.cells.begin(); cell < W.cells.end(); cell++)
    {
  	    // foreach cell go through it's particles... 
	    for (std::vector<Particle>::iterator i = cell->particles.begin(); i < cell->particles.end(); i++)
        {
            // compare Cell Numbers of (maybe) moved particle i
            int checkCell = W.getCellNumber(i);
            // 	..at first calc new position in every dimension
            for (unsigned int d=0; d<DIM; d++)
		    {
                // computing new location of the particle i if it's leaving the world, elsewise just call handle_borders (-lc version) in the end
	  	        i->x[d] += W.delta_t*i->v[d] + (.5*i->F[d]*sqr(W.delta_t)) / i->m;

                // DEBUG:
//                std::cout << "Cell[" << cell-W.cells.begin() << "]"
//                          << ".particles["  <<  i-cell->particles.begin() << "]"
//                          << ".x[" << d << "]=" << i->x[d] << std::endl;
 
                // save last force...
		        i->F_old[d] = i->F[d];
                // ... and don't forget to set the actual force to zero
                i->F[d] = 0;
		    }

            // then check if particle left its cell and handle moving issues (respective border issues)
            if (W.getCellNumber(i) != checkCell)
            {
                // if the flag is checked, push the particle in the last round into it's new position
                doIt = false;
                innerWorld = true;
                // check it in every dimension
                for (unsigned int d = 0; d<DIM; d++)
                {
                    // periodic - position = position % worldlength
                    if (i->x[d] > W.length[d] && W.upper_border[d] == W.periodic)
                    {
                        // new position is at the beginning of world plus the overhead that x left the world
                        i->x[d] = fmod(i->x[d], W.length[d]);
                        doIt = true;
                        innerWorld = false;
                    }
                    // periodic - position = position % worldlength
                    else if (i->x[d] < 0 && W.lower_border[d] == W.periodic)
                    {
                        // new position is end of world minus overhead that x left the world
                        i->x[d] =  W.length[d] - fabs(fmod(i->x[d], W.length[d]));
                        doIt = true;
                        innerWorld = false;
                    }
                    // leaving - it just bumps out
                    else if (i->x[d] > W.length[d] && W.upper_border[d] == W.leaving)
                    {
                        innerWorld = false;
                        // DEBUG:
                        std::cout << "New position (oben raus Wegvomfenster): " << std::endl;
                        // regardless in which dimension, just erase
                        i = cell->particles.erase(i);
                        i--;
                        W.nParticles--;
                        // don't forget to set dimension to DIM  or you will handle another particle in the wrong dimension
                        d=DIM;
                        break;
                    }
                    // leaving - it just bumps out
                    else if (i->x[d] < 0  && W.lower_border[d] == W.leaving)
                    {
                        innerWorld = false;
                        // DEBUG:
                        std::cout << "New position (unten raus Wegvomfenster): " << std::endl;
                        // regardless in which dimension, just erase
                        i = cell->particles.erase(i);
                        i--;
                        W.nParticles--;
                        d=DIM;
                        break;
                    }
                    // the particle just changes the cell
                    //else

                    if(innerWorld)
                    {
                        int tmp = W.getCellNumber(i);
                        int tmp2 = W.cells.size();
                        //W.cells[W.getCellNumber(i)].particles.push_back(cell->particles[i-cell->particles.begin()]);
                        // push the moved cell temporary to worlds particle list and reinsert it later
                        W.particles.push_back(*i);
                        i = cell->particles.erase(i);
                        i--;
                        d=DIM;
                        break;
                    }
                    else
                        std::cout << " PROBLEM PROBLEM PROBLEM PROBLEM PROBLEM PROBLEM PROBLEM PROBLEM " << std::endl;

                    // in the last round, push the particles when it left periodically in whatever dimension
                    if (d==DIM-1 && doIt)
                    {
                        // DEBUG:
//                            std::cout << "New position (oben/unten raus Periodisch): " << i->x[d] << std::endl;

                        // Add particle i to its corresponding cell
                        W.cells[W.getCellNumber(i)].particles.push_back(*i);
                        i = cell->particles.erase(i);
                        i--;
                        break;
                    }
                }
            }
        }
    }
    // and now add the particles again to their belonging cells
    for (std::vector<Particle>::iterator i = W.particles.begin(); i < W.particles.end(); i++)
    {
        //W.cells[W.getCellNumber(i)].particles.push_back(cell->particles[i-cell->particles.begin()]);
//        W.cells[W.getCellNumber(i)].particles.push_back(W.particles[i-W.particles.begin()]);
        if (W.getCellNumber(i)<0) {
        } else {
        W.cells[W.getCellNumber(i)].particles.push_back(*i);
        i = W.particles.erase(i);
        i--;
        }
    }
}
// vim:set et sts=4 ts=4 sw=4 ai ci cin:
