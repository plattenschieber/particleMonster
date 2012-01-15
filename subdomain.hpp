#ifndef SUBDOMAIN_HPP
#define SUBDOMAIN_HPP

#include "defines.hpp"


/**
 * @brief Particle data
 *
 * this class contains the particle data
 */
class SubDomain
{
public:
    /// edge length of domain
    real L[DIM];
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
    int N_p[DIM];

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
};

#endif // SUBDOMAIN_HPP
// vim:set et sts=4 ts=4 sw=4 ai ci cin:
