#include "velocityverlet.hpp"
#include <math.h> 
VelocityVerlet::VelocityVerlet(World& _W, Potential& _Pot, ParObserverXYZ& _O) : TimeDiscretization(_W,_Pot,_O)
{
    // empty constructor(CBR) - not really, inheritance of TimeDiscretization sets all variables by his initializer list
}

VelocityVerlet::VelocityVerlet(World& _W, Potential* _Pot, ParObserverXYZ& _O) : TimeDiscretization(_W,(*_Pot),_O)
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
    // increase time
    W.t += delta_t;
    //if (W.s.myrank == 0)
        std::cout << "STEP " << W.t << std::endl;
    // count steps
    W.step++;
    // first of all, we have to update our positions by means of their actual position, pace and force
    updateX();
    // then we update their new force, based on the potential and the new world situation
    compF();
    // now we can compute their new pace
    updateV();
    // communicate energy
    MPI_Allreduce(&W.e_kin, &W.e_kin, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&W.e_pot, &W.e_pot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&W.e_tot, &W.e_tot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    // compute new energy average
    updateAverage();
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
    // velocity scaling
    real beta = 1.0;

    // roll over every particle i...
    for (std::list<Particle>::iterator i = W.particles.begin(); i != W.particles.end(); i++)
        // ...and every of it's dimensions
        for (unsigned int d=0; d<DIM; d++)
        {
            // compute new velocity in dimension d
            i->v[d] += .5*(i->F_old[d] + i->F[d])*W.delta_t/i->m;
            // add now the pro rata e_kin
            W.e_kin += .5*i->m*sqr(i->v[d]);

            // VELOTCITY SCALING
            if (W.isThermoStartTemp && (W.step % W.T_Step == 0) && (fabs(W.T_D - W.T) > 10e-6) )
            {
                beta = sqrt(W.T_D * (W.nParticles-1) / (48*W.e_kin));
                std::cout << "VELOTCITY SCALING" << std::endl;
                // multiply velocity by beta
                i->v[d] *= beta;
                // and scale kinetc energy
                W.e_kin *= beta;
            }
        }
    // compute total energy
    W.e_tot = W.e_kin + W.e_pot;
    // set new temperature
    W.T *= sqr(beta);
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
            if ( (W.upper_border[d] == W.leaving) & (i->x[d]>0) & (i->x[d] > W.worldLength[d]) )
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

void VelocityVerlet::updateAverage()
{
    // add new energy to our average list
    W.ekin_list.push_back (W.e_kin);
    W.epot_list.push_back (W.e_pot);
    // if there are more than 100 measurements, pop the oldest (the first) and update sums
    if (W.step > 100)
    {
        W.ekin_sum -= W.ekin_list.front ();
        W.epot_sum -= W.epot_list.front ();
        W.ekin_list.pop_front ();
        W.epot_list.pop_front ();
    }
    W.ekin_sum += W.e_kin;
    W.epot_sum += W.e_pot;
    // Divide by number of elements in list or at most 100
    W.ekin_avg = W.ekin_sum / ((W.step<100)?W.step:100);
    W.epot_avg = W.epot_sum / ((W.step<100)?W.step:100);

}


// vim:set et sts=4 ts=4 sw=4 ai ci cin:
