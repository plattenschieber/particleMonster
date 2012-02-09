#ifndef _SUBDOMAIN_HPP
#define _SUBDOMAIN_HPP

#include "defines.hpp"
#include "worldlc.hpp"
#include "cell.hpp"
#include "defines.hpp"
#include "particle.hpp"
#include <mpi.h>
#include <vector>
#include <list>
#include <string>
#include <iostream>
#include <fstream>
#include <map>

/**
 * @brief Particle data
 *
 * this class contains the particle data
 */
class SubDomain : public WorldLC {
public:
    /// standard ctor
    SubDomain();
    /* @brief Overwrite existing World::readParameter to handle new files
     *
     * parameter file example TODO: rewrite code example
     * \code
     * \endcode
     *
     * @param filename filename of the parameter file
     */
    void readParameter (const std::string &filename);
    // TODO: Add some comment here
    void readParticles (const std::string &filename);

    int getCellNumber (const Particle &p);

    void sendReceive (int lower_proc, int *lower_ic_start,  int *lower_ic_stop, int *lower_ic_startreceive, int* lower_ic_stopreceive,
                      int upper_proc, int *upper_ic_start,  int *upper_ic_stop, int *upper_ic_startreceive, int *upper_ic_stopreceive);

    void deleteBorderParticles ();
    void constructParticle (MPI::Datatype& MPI_Particle);

    void setCommunication ( int d,
                            int *lower_ic_start, int *lower_ic_stop, int *lower_ic_startreceive, int *lower_ic_stopreceive,
                            int *upper_ic_start, int *upper_ic_stop, int *upper_ic_startreceive, int *upper_ic_stopreceive);


    void communicate (bool isForward);

    /// edge length of domain
    real domainLength[DIM];
    /// number of cells in the whole domain
    int N_c[DIM];

    // further parameter for parallelization
    /// process rank of local process
    int myrank;
    //// number of started processes
    int numprocs;
    /// process position in grid
    int ip[DIM];
    /// number of subdomains
    int nProcs[DIM];

    /// process number of neighbourprocess
    int ip_lower[DIM];
    int ip_upper[DIM];
    /// width of bordure, correspond to to first inner local
    int ic_start[DIM];
    int ic_stop[DIM];
    /// number of cells in subdomain with bordure
    int ic_number[DIM];
    /// cell length
    real cellh[DIM];
    /// global index of first cell in sub field
    int ic_lower_global[DIM];
    /// global index of last cell in sub field
    int ic_upper_global[DIM];
    /// a particle for MPI
    MPI::Datatype MPI_Particle;
};

/**
 * @brief a ostream operator for the SubDomain class
 *
 * @param os stream object
 * @param S the SubDomain
 *
 * @return resulting stream object
 */
std::ostream& operator << (std::ostream& os, SubDomain& S);

#endif // _SUBDOMAIN_HPP
// vim:set et sts=4 ts=4 sw=4 ai ci cin:
