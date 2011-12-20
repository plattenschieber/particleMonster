#include "velocityverlet.hpp"
#include <math.h> 
VelocityVerlet::VelocityVerlet(World& _W, Potential& _Pot, ObserverXYZ& _O) : TimeDiscretization(_W,_Pot,_O)
{
    // empty constructor(CBR) - not really, inheritance of TimeDiscretization sets all variables by his initializer list
}

VelocityVerlet::VelocityVerlet(World& _W, Potential* _Pot, ObserverXYZ& _O) : TimeDiscretization(_W,(*_Pot),_O)
{
    // empty constructor(CBV) - not really, inheritance of TimeDiscretization sets all variables by his initializer list
}

// runs the simulation
void VelocityVerlet::simulate()
{
    // while simulation end time not reached
    while (W.t < W.tEnd)
        // it's just one step for a man, but a big step for the world...
        timestep(W.deltaT);
}

// a timestep is a timestep is a timestep. It doesn't matter, when it happens. Our model stays consistent!
void VelocityVerlet::timestep(real deltaT)
{
    // first of all, we have to update our positions by means of their actual position, pace and force
    updateX();
    // then we update their new force, based on the potential and the new world situation
    compF();
    // now we can compute their new pace
    updateV();
    
    // if they left the world, no other treatment of the particles is neaded
    //handleBorders();
    // increase time
    W.t += deltaT;
    // notify observer
    O.notify();
}

void VelocityVerlet::compF()
{
    // check the distance and throw out all that's more far away than rcut
    real dist = 0.0, rcut = 2.5;
    // there is no ePot in the beginning
    W.ePot = 0.0;
    // we compute the ePot for each pair of particles and add it to the worlds' ePot...
    for (std::list<Particle>::iterator i = W.particles.begin(); i != W.particles.end(); i++)
        // ...except of the computation with itself (i!=j)
        for (std::list<Particle>::iterator j = W.particles.begin(); j != i; j++)
        {
            // don't forget to reset the distance
            dist = 0.0;
            // check the distance
    		for (int d=0; d<DIM; d++)
		    dist += sqr(j->x[d]-i->x[d]);	
		    // only particles which are closer than rcut
		    if(dist <= rcut) 
		         // computes the force between particle i and j and add it to our potential
             W.ePot += Pot.force(*i, *j);
		}
}

void VelocityVerlet::updateV()
{
    // there is no eKin in the beginning
    W.eKin = 0.0;
	// roll over every particle i...
    for (std::list<Particle>::iterator i = W.particles.begin(); i != W.particles.end(); i++)
        // ...and every of it's dimensions
		for (unsigned int d=0; d<DIM; d++)
        {
            // compute new velocity in dimension d
            i->v[d] += .5*(i->F_old[d] + i->F[d])*W.deltaT/i->m;
            // add now the pro rata eKin
            W.eKin += .5*i->m*sqr(i->v[d]);
        }
}

void VelocityVerlet::updateX()
{
    // roll over every particle...
    for (std::list<Particle>::iterator i = W.particles.begin(); i != W.particles.end(); i++)
        // ...and every of it's dimensions
		for (unsigned int d=0; d<DIM; d++)
		{
		    // computing new location of the particle i
            i->x[d] += W.deltaT*i->v[d] + (.5*i->F[d]*sqr(W.deltaT)) / i->m;
		    // save last force...
			i->F_old[d] = i->F[d];
		    // ... and don't forget to set the actual force to zero
			i->F[d] = 0;
		}
}

void VelocityVerlet::handleBorders()
{
    // roll over every particle...
    for (std::list<Particle>::iterator i = W.particles.begin(); i != W.particles.end(); i++)
        // ...and every of it's dimensions
		for (unsigned int d=0; d<DIM; d++)
		{
			// is leaving AND above zero AND outer space
                    if ( (W.upperBorder[d] == W.leaving) & (i->x[d]>0) & (i->x[d] > W.length[d]) )
		        {
                      W.particles.erase(i);
                      i--;
                      break;
		        }
            		// same here, except of handling the particles under zero
                    else if( (W.lowerBorder[d] == W.leaving) && (i->x[d]<0) )
            		{
                        W.particles.erase(i);
                        i--;
                        break;
           		}
            		// all other particles are at least on the point zero. This point is actually included ;) 
       		 }
}

// vim:set et sts=4 ts=4 sw=4 ai ci cin:
