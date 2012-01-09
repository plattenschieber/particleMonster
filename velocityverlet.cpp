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
    while (W.t < W.t_end)
        // it's just one step for a man, but a big step for the world...
        timestep(W.delta_t);
}

// a timestep is a timestep is a timestep. It doesn't matter, when it happens. Our model stays consistent!
void VelocityVerlet::timestep(real delta_t)
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
    W.t += delta_t;
<<<<<<< HEAD

    std::cout << "STEP " << W.t << std::endl;

=======
>>>>>>> parent of 0a679d7... Revert 9c405afdfa97c3bcde5d29dbfef7d7b9ff9c8dd9^..HEAD
    // notify observer
    O.notify();
}

void VelocityVerlet::compF()
{
    // check the distance and throw out all that's more far away than rcut
    real dist = 0.0, rcut = 2.5;
    // there is no e_pot in the beginning
    W.e_pot = 0.0;
    // we compute the e_pot for each pair of particles and add it to the worlds' e_pot...
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
                W.e_pot += Pot.force(*i, *j);
        }
}

void VelocityVerlet::updateV()
{
    // there is no e_kin in the beginning
    W.e_kin = 0.0;
    // roll over every particle i...
    for (std::list<Particle>::iterator i = W.particles.begin(); i != W.particles.end(); i++)
        // ...and every of it's dimensions
        for (unsigned int d=0; d<DIM; d++)
        {
            // compute new velocity in dimension d
            i->v[d] += .5*(i->F_old[d] + i->F[d])*W.delta_t/i->m;
<<<<<<< HEAD
            // if we want to check the temperatur regulary and there is a Thermo option
            if (fmod(W.t,W.thermo_step_interval) == 0 && W.isThermoStartTemp)
                // multiply velocity by beta
                i->v[d] *= W.calcBeta();
=======
>>>>>>> parent of 0a679d7... Revert 9c405afdfa97c3bcde5d29dbfef7d7b9ff9c8dd9^..HEAD
            // add now the pro rata e_kin
            W.e_kin += .5*i->m*sqr(i->v[d]);
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
            i->x[d] += W.delta_t*i->v[d] + (.5*i->F[d]*sqr(W.delta_t)) / i->m;
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
            if ( (W.upper_border[d] == W.leaving) & (i->x[d]>0) & (i->x[d] > W.length[d]) )
            {
                W.particles.erase(i);
                i--;
                break;
            }
            // same here, except of handling the particles under zero
            else if( (W.lower_border[d] == W.leaving) && (i->x[d]<0) )
            {
                W.particles.erase(i);
                i--;
                break;
            }
            // all other particles are at least on the point zero. This point is actually included ;)
        }
}

// vim:set et sts=4 ts=4 sw=4 ai ci cin:
