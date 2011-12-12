// in addiction to our DIM we will span up a "for"-tree
#include "velocityverletlc.hpp"
#include <math.h> 
#include "defines.hpp"



// another macro to expand - this time a for loop based on dimension DIM
#

VelocityVerletLC::VelocityVerletLC(WorldLC& _W, LJPotential& _Pot, ObserverXYZ& _O) : VelocityVerlet(_W,_Pot,_O)
{
    // initialize your own World, otherwise implicit cast to World will force us to explicit cast the World every time we use it, to WorldLC
    W = _W;
    Pot = _Pot;
}
VelocityVerletLC::VelocityVerletLC(WorldLC& _W, LJPotential* _Pot, ObserverXYZ& _O) : VelocityVerlet(_W,(*_Pot),_O)
{
    // initialize your own World, otherwise implicit cast to World will force us to explicit cast the World every time we use it, to WorldLC
    W = _W;
    Pot =(* _Pot);
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
		  // don't forget to reset the distance
		  dist = 0.0;
		  for (int d=0; d<DIM; d++)
		     // accumulate distance between particle and neighbour cell. 
		     dist += sqr(W.cell_length[d]*nbCell[d] - i->x[d]);	
		     // periodic , unknown , leaving  Handle borders more specific
		     if (nbCell[d]<0 && W.lower_border[d]==W.periodic) nbCell[d]=W.cell_N[d]; 
		     else if (nbCell[d]>W.cell_N[d] && W.upper_border[d]==W.periodic) nbCell[d]=0; 
		  // If neighbour cell is too far away, no calc of inner (more far away) particles are needed! 
		  if (dist <= W.cell_r_cut)
	             for (std::vector<Particle>::iterator j = W.cells[J(jCell,W.cell_N)].particles.begin(); j < W.cells[J(jCell,W.cell_N)].particles.end();)
                        // ...except of the computation with itself (i!=j) 
			if(i!=j)
			    // don't forget to reset the distance
			    dist = 0.0;
			    // check the distance
			    for (int d=0; d<DIM; d++) dist += sqr(j->x[d] - i->x[d]);	
				// only particles which are closer than rcut
				if(dist <= W.cell_r_cut ) 
				     // computes the force between particle i and j and add it to our potential
				     W.e_pot += Pot.force(*i, *j, dist, W.epsilon, W.sigma);
    {
        for (jCell[1]=0; jCell[1]<W.cell_N[1]; jCell[1]++)
        {
        	for (jCell[2]=0; jCell[2]<W.cell_N[2]; jCell[2]++)
            {
                    for (nbCell[0]=jCell[0]-1; nbCell[0]<=jCell[0]+1; nbCell[0]++)
                    {
                        for (nbCell[1]=jCell[1]-1; nbCell[1]<=jCell[1]+1; nbCell[1]++)
                        {
                           for (nbCell[2]=jCell[2]-1; nbCell[2]<=jCell[2]+1; nbCell[2]++)
                           {
                              {
                              }
                              { 
                                {
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
	// roll over every cell	
   	for (std::vector<Cell>::iterator cell =  W.cells.begin(); cell < W.cells.end(); cell++)
    {
  	    // foreach cell go through it's particles... 
	    for (std::vector<Particle>::iterator i = cell->particles.begin(); i < cell->particles.end(); i++)
        {
            // ...and over every dimension of particle i
		    for (unsigned int d=0; d<DIM; d++)
		    {
                // DEBUG:
                std::cout << "Aktuelle Position von Partikel " <<  i-cell->particles.begin() << " << i->x[" << d << "]=" << i->x[d] << std::endl;
                // computing new location of the particle i if it's leaving the world, elsewise just call handle_borders (-lc version) in the end
	  	        i->x[d] += W.delta_t*i->v[d] + (.5*i->F[d]*sqr(W.delta_t)) / i->m;
                std::cout << "Neue Position von Partikel " <<  i-cell->particles.begin() << " << i->x[" << d << "]=" << i->x[d] << std::endl;
		        // periodic - position = position % worldlength
		        if (i->x[d] > W.length[d] && W.upper_border[d] == W.periodic)
                {
                    i->x[d] = fmod(i->x[d], W.length[d]);
                    // DEBUG:
                    std::cout << "Neue Position (oben raus Periodisch): " << i->x[d] << std::endl;

                }
		        if (i->x[d] < 0 && W.lower_border[d] == W.periodic)
                {
		            i->x[d] =  W.length[d] - fabs(fmod(i->x[d], W.length[d]));
		            // DEBUG:
                    std::cout << "Neue Position (unten raus Periodisch): " << i->x[d] << std::endl;
                }
		        // leaving - it just bumps out
		        if (i->x[d] > W.length[d] && W.upper_border[d] == W.leaving) 
		        {
		            // DEBUG:
                    std::cout << "Neue Position (oben raus Wegvomfenster): " << std::endl;
		            cell->particles.erase(i); d=DIM; break; 
                }
		        if (i->x[d] < 0  && W.lower_border[d] == W.leaving)
		        { 
		            // DEBUG:
                    std::cout << "Neue Position (unten raus Wegvomfenster): " << std::endl;
		            cell->particles.erase(i); d=DIM; break; 
		        }
                // save last force...
		        i->F_old[d] = i->F[d];
                // ... and don't forget to set the actual force to zero
	    	    i->F[d] = 0;
		    }
        }
    }
}
// put the particles into their right places
void VelocityVerletLC::handle_borders()
{
    int jCell[DIM], nbCell[DIM];
    // roll over every cell
    for (jCell[0]=0; jCell[0]<W.cell_N[0]; jCell[0]++)
    {
    	for (jCell[1]=0; jCell[1]<W.cell_N[1]; jCell[1]++)
        {
            for (jCell[2]=0; jCell[2]<W.cell_N[2]; jCell[2]++)
	        {
	            for (std::vector<Particle>::iterator i = W.cells[J(jCell,nbCell)].particles.begin(); i < W.cells[J(jCell,nbCell)].particles.end(); i++)   
                {
                    for (int d=0; d<DIM; d++)
			            nbCell[d] = (int)floor(i->x[d] * W.cell_N[d] / W.cell_length[d]);
		            if((jCell[0]!=nbCell[0] || jCell[1]!=nbCell[1] || jCell[2]!=nbCell[2]))
            			W.cells[J(nbCell,W.cell_N)]; 
                }
            }
        }
    }   
}
// vim:set et sts=4 ts=4 sw=4 ai ci cin:
