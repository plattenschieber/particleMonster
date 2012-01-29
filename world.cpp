#include "world.hpp"



// intializing via intializer list - order is set by definings in .hpp;
World::World() : name("unknown"),t(0),delta_t(0.0),t_end(0.0),e_kin(0.0),e_pot(0.0),e_tot(0.0),e_avg(0.0),
    sigma(1.0),epsilon(1.0), nParticles(0),
    thermo_start_temp(0.0),isThermoStartTemp(false),
    thermo_step_interval(0.0),thermo_target_temp(0.0)
{
    //setting map mapOptions for switchcase
    mapOptions["name"] = NAME;
    mapOptions["delta_t"] = DELTA_T;
    mapOptions["t_end"] = T_END;
    mapOptions["length"] = LENGTH;
    mapOptions["upper_border"] = UPPERBORDER;
    mapOptions["lower_border"] = LOWERBORDER;
    mapOptions["sigma"] = SIGMA;
    mapOptions["epsilon"] = EPSILON;
    mapOptions["set_start_temperature"] = STARTTEMP;
    mapOptions["thermostat_step_interval"] = STEPINTERVAL;
    mapOptions["thermostat_target_temperature"] = TARGETTEMP;
    mapOptions["random_seed"] = RANDOMSEED;
}

void World::readParameter(const std::string &filename)
{

    // create input filestream
    std::ifstream parfile(filename.c_str());
    // check if file is open
    if (!parfile.is_open())
        throw std::runtime_error("readParameter(): Can't open file '" + filename + "' for reading.");

    // helper strings
    std::string line, option, tmp;

    // read file till eof
    while (parfile.good())
    {
        // read line from file
        getline(parfile,line);
        // create a string stream
        std::stringstream strstr;
        // put line into string stream
        strstr << line;
        // read option and value from stringstream
        strstr >> option;
        // push next read value, with internal converter of string stream, into the propper place
        switch(mapOptions[option])
        {
        case DELTA_T:
            strstr >> delta_t;
            break;
        case T_END:
            strstr >> t_end;
            break;
        case NAME:
            strstr >> name;
            break;
        case EPSILON:
            strstr >> epsilon;
            break;
        case SIGMA:
            strstr >> sigma;
            break;
        case LENGTH:
            for (int i=0; i<DIM; i++) {
                strstr >> length[i];
            }
            break;
        case UPPERBORDER:
            strstr >> tmp;
            for (int i=0; i<DIM; i++)
            {
                if (tmp == "leaving") upper_border[i] = leaving;
                else if (tmp == "periodic") upper_border[i] = periodic;
                else upper_border[i] = unknown;
            }
            break;
        case LOWERBORDER:
            strstr >> tmp;
            for (int i=0; i<DIM; i++)
            {
                if (tmp == "leaving") lower_border[i] = leaving;
                else if (tmp == "periodic") lower_border[i] = periodic;
                else lower_border[i] = unknown;
            }
            break;

        case STARTTEMP:
            strstr >> thermo_start_temp;
            isThermoStartTemp = true;
            break;
        case STEPINTERVAL:
            strstr >> thermo_step_interval;
            break;
        case TARGETTEMP:
            strstr >> thermo_target_temp;
            break;
        case RANDOMSEED:
            double tmp2;
            strstr >> tmp2;
            if (tmp2<1)
                srand(time(NULL));
            else
                srand(tmp2);
            break;
            // handle unknown options
        default:
            std::cout << "'" << option << "' is an invalid option." << std::endl;
            break;
        }
    }
    // close file
    parfile.close();
}

void World::readParticles(const std::string &filename)
{
    // create input filestream
    std::ifstream parfile(filename.c_str());
    // check if file is open
    if (!parfile.is_open())
        throw std::runtime_error("readParticles(): Can't open file '" + filename + "' for reading.");

    // helper strings
    std::string line;
    Particle tmpparticle;

    // read file till eof
    while (parfile.good())
    {
        // read line from file
        getline(parfile,line);
        // break if there is a blank line, e.g. a newline at the eof
        if (line=="")
            break;
        // create a string stream
        std::stringstream strstr;
        // push line into string stream
        strstr << line;
        // grab ID:
        strstr >> tmpparticle.ID;
        // save the particles mass
        strstr >> tmpparticle.m;
        // push dim read values into tmpparticles position and check wether particle is inside the world, or not
        for(int d=0; d<DIM; d++)
        {
            strstr >> tmpparticle.x[d];
            if (tmpparticle.x[d] < 0 || tmpparticle.x[d] > length[d])
            {
                std::cerr << "<------- FAILURE ------->" << std::endl;
                std::cerr << "There is a particle outside the world. Please check input file." << std::endl;
                exit (EXIT_FAILURE);
            }

        }

        // push dim read values into tmpparticles velocity
        for(int d=0; d<DIM; d++)
            strstr >> tmpparticle.v[d];
        // set start velocity of particle by means of Maxwell Boltzmann
        if (isThermoStartTemp)
        {
            // number crunchers for Maxwell-Boltzmann
            real s, r, u[DIM], tmp;
            do
            {
                // reset
                tmp = 0.0;
                for (int d=0; d<DIM; d++)
                {
                    // get out a number between 0.0 and 1.0 inclusively
                    u[d] = (double(rand()) / RAND_MAX);
                    // scaling
                    u[d] = 2.0*u[d]-1.0;
                    tmp += sqr(u[d]);
                }
                s = tmp;

                // are we inside the R^DIM ball? sin^2 + cos^2 = 1
            } while (s>1.0);
            // i don't know how it's working exactly
            r = -2.0*log(s)/s;
            // set new velocity according to Maxwell-Boltzmann
            for (int d=0; d<DIM; d++)
                tmpparticle.v[d] = u[d]*sqrt(r);
            tmpparticle.v[2] = 0;
        }

        // assure integrity of the force of our tmpparticle
        for(int d=0; d<DIM; d++)
            tmpparticle.F[d] = tmpparticle.F_old[d] = 0.0;

        // calc according cell
        int cell[DIM];
        // is particle in our Subdomain
        bool isInSubdomain = true;
        // check if particle doesn't belong to our subdomain
        for (int d=0; d<DIM; d++)
        {
            // calc cell number and offset
            cell[d] = (int)floor( tmpparticle.x[d]/s.cellh[d] ) - s.ic_lower_global[d] + s.ic_start[d];
            if( cell[d] < s.ic_start[d] || cell[d] >= s.ic_stop[d] )
                isInSubdomain = false;
        }

        if (isInSubdomain)
        {
            // add the new particle to our worlds' particles
            particles.push_back(tmpparticle);
            // update particle size
            nParticles++;
        }
    }
    // close file
    parfile.close();
}


real World::calcBeta()
{
    return sqrt(thermo_target_temp * (nParticles-1) / (48*e_kin));
}


std::ostream& operator << (std::ostream& os, World& W) {
    os << "t=" << W.t << " delta_t=" << W.delta_t << " t_end=" << W.t_end
       << " Number of Particles=" << W.particles.size();
    for (std::list<Particle>::iterator i = W.particles.begin(); i != W.particles.end(); i++)
    {
        os << "particle[" << i->ID << "] = (";
        for (unsigned int d=0; d<DIM; d++)
            os << i->x[d] << ", ";
        os << ")";
    }
    return os;
}
// vim:set et sts=4 ts=4 sw=4 ai ci cin:
