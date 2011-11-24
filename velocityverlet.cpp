#include "velocityverlet.hpp"
#include <math.h> 
VelocityVerlet::VelocityVerlet(World& _W, Potential& _Pot, Observer& _O) : TimeDiscretization(_W,_Pot,_O)
{
    // empty constructor(CBR) - not really, inheritance of TimeDiscretization sets all variables by his initializer list
}

VelocityVerlet::VelocityVerlet(World& _W, Potential* _Pot, Observer& _O) : TimeDiscretization(_W,(*_Pot),_O)
{
    // empty constructor(CBV) - not really, inheritance of TimeDiscretization sets all variables by his initializer list
}

// runs the simulation
void VelocityVerlet::simulate()
{
    // while simulation end time not reached
    while (W.t < W.t_end)
        // it's just one step for a man, but a big step for the world...
        timestep(W.delta_t);
}

// a timestep is a timestep is a timestep. It doesn't matter, when it happens. Our model stays consistent!
void VelocityVerlet::timestep(real delta_t)
{
    // first of all, we have to update our positions by means of their actual position, pace and force
	update_X(); 
    // then we update their new force, based on the potential and the new world situation
	comp_F();
    // now we can compute their new pace
	update_V();
    
    // if they left the world, no other treatment of the particles is neaded
    handle_borders();
    // increase time
    W.t += delta_t;
    // notify observer
    O.notify();
}

void VelocityVerlet::comp_F()
{
	// check the distance and throw out all that's more far away than rcut
	real tmp = 0.0, rcut = 2.5;
    // there is no e_pot in the beginning
	W.e_pot = 0.0;
    // we compute the e_pot for each pair of particles and add it to the worlds' e_pot...
	for (std::vector<Particle>::iterator i = W.particles.begin(); i < W.particles.end(); i++)
        // ...except of the computation with itself (i!=j)
		for (std::vector<Particle>::iterator j = W.particles.begin(); j < i; j++) // changed from W.particles.size() && i!=j and removed the *.5 in the inner for-loop // unsigned int j=0; j<i; j++
		{
			tmp = 0.0;
			// check the distance
	    		for (int d=0; d<DIM; d++)
				tmp += sqr(j->x[d]-i->x[d]);	
			// only particles which are closer than rcut
			if(tmp<=rcut) 
				// computes the force between particle i and j and add it to our potential
				 W.e_pot += Pot.force(*i, *j);
		}
}

void VelocityVerlet::update_V()
{
    // there is no e_kin in the beginning
	W.e_kin = 0.0;
	// roll over every particle i...
    for (std::vector<Particle>::iterator i = W.particles.begin(); i < W.particles.end(); i++)
        // ...and every of it's dimensions
		for (unsigned int d=0; d<DIM; d++)
        {
            // compute new velocity in dimension d
			i->v[d] += .5*(i->F_old[d] + i->F[d])*W.delta_t/i->m;
            // add now the pro rata e_kin 
            W.e_kin += .5*i->m*sqr(i->v[d]);
        }
}

void VelocityVerlet::update_X()
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

void VelocityVerlet::handle_borders()
{
    // roll over every particle...
    for (std::vector<Particle>::iterator i = W.particles.begin(); i < W.particles.end(); i++)
        // ...and every ojf it's dimensions
		for (unsigned int d=0; d<DIM; d++)
		{
            // if particle is above zero, there are no borders and it's further than our worlds length in this axis, then pop it
               // DEBUG  std::string tmp; 
                // DEBUG if(i->x[d] > W.length[d])
                   // DEBUG  tmp ="hans"; else tmp ="peter";
               // DEBUG  std::cout << tmp << W.length[d] << "und" << i->x[d] << "welt" << W.particles.size() <<std::endl;
            if ( (W.upper_border[d] == W.leaving) &&
                (i->x[d]>0) &&
                (i->x[d] > W.length[d]) )
            {
                W.particles.erase(i);
                break;
            }
            // same here, except of handling the particles under zero
            else if( (W.lower_border[d] == W.leaving) && (i->x[d]<0) && (i->x[d] < -1* W.length[d]) )
            {
                W.particles.erase(i);
                break;
            }
            // all other particles are at least on the point zero. This point is actually included ;) 
        }
}

// vim:set et sts=4 ts=4 sw=4 ai ci cin:
