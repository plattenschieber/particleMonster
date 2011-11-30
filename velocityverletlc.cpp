#ifndef _WORLDLC_HPP
#define _WORLDLC_HPP

// index J of cell will be computed by our cell indices jCell
#if DIM==2
#define J(jCell,cell_N) ((jCell)[0] + (cell_N)[0]*(jCell)[1]) 
#elif DIM==3
#define J(jCell,cell_N) ((jCell)[0] + (cell_N)[0]*((jCell[1] + (cell_N)[1]*(jCell)[2])))
#endif

// in addiction to our DIM we will span up a "for"-tree
#if DIM==2
#define foreachdim(iCell
#include "velocityverletlc.hpp"
#include <math.h> 
// TODO: Is this already right, with the initializer list?
VelocityVerletLC::VelocityVerletLC(WorldLC& _W, Potential& _Pot, Observer& _O) : VelocityVerlet(_W,_Pot,_O)
{
    // initialize your own World, otherwise implicit cast to World will force us to explicit cast the World every time we use it, to WorldLC
    W = _W;
}
// TODO: same here!!
VelocityVerletLC::VelocityVerletLC(WorldLC& _W, Potential* _Pot, Observer& _O) : VelocityVerlet(_W,(*_Pot),_O)
{
    // initialize your own World, otherwise implicit cast to World will force us to explicit cast the World every time we use it, to WorldLC
    W = _W;
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
    for (std::vector<Cell>::iterator cell = W.cells.begin(); cell < W.cells.end(); cell++)
        // we compute the e_pot for each pair of particles in it's cell including the neighbour cells and add it to the worlds' e_pot...
	// roll over every particle i in actual cell
        for (std::vector<Particle>::iterator i = cell->particles.begin(); i < cell->particles.end(); i++)
	{
	    // set all F_i to zero
	    for (int d=0; d<DIM; d++) i->F[d]=0;
	    // roll over every neighbour cell
	    for (std::vector<Cell>::iterator neighbour = cell+1; neighbour += 2)    
                // ...except of the computation with itself (i!=j)
	        for (std::vector<Particle>::iterator j = neighbour->particles.begin(); j < neighbour->particles.end(); j++) 
	        {
			if(i!=j)
			{
			    // don't forget to reset the distance
			    dist = 0.0;
			    // check the distance
			    for (int d=0; d<DIM; d++) dist += sqr(j->x[d] - i->x[d]);	
				// only particles which are closer than rcut
				if(dist <= rcut) 
				     // computes the force between particle i and j and add it to our potential
				     W.e_pot += Pot.force(*i, *j);
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
		// foreach cell go through it's particles... 
		for (std::vector<Particle>::iterator i = cell->particles.begin(); i < cell->particles.end(); i++)
	    	// ...and over every dimension of particle i
	    		for (unsigned int d=0; d<DIM; d++)
            		{
				// compute new velocity in dimension d
				i->v[d] += .5*(i->F_old[d] + i->F[d])*W.delta_t/i->m;
                		// add now the pro rata e_kin 
                		W.e_kin += .5*i->m*sqr(i->v[d]);
            		}


}

void VelocityVerletLC::update_X()
{
    // roll over every particle...
    for (std::vector<Particle>::iterator i = W.particles.begin(); i < W.particles.end(); i++)
        // ...and every of it's dimensions
		for (unsigned int d=0; d<DIM; d++)
		{
            // computing new location of the particle i
			i->x[d] += W.delta_t*i->v[d] + (.5*i->F[d]*sqr(W.delta_t)) / i->m;
            // save last force...
			i->F_old[d] = i->F[d];
            // ... and don't forget to set the actual force to zero
			i->F[d] = 0;
		}
}

// vim:set et sts=4 ts=4 sw=4 ai ci cin:
