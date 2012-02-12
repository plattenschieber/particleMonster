/** \mainpage MD Praktikum 2011
 *
 * \section blatt1_sec Übungsblatt 1
 *
 * Simulation von Planeten- und Kometenbahnen.
 *
 * \subsection algorithms_subsec Implementierte Verfahren und deren Klassen
 *
 * - GravityPotential, das (skalierte) Gravitationspotential \f$U(r_{ij}) = -m_im_j/r_{ij}\f$
 * - LJPotential, ein Pseudopotential der klassischen Moleküldynamik \f$U(r_{ij}) = 4 \cdot \epsilon (\frac{\sigma}{r_{ij}})^6 \cdot    ( (\frac{\sigma}{r_{ij}})^6 - 1 ) \f$
 * - VelocityVerlet mit \f$O(N^2)\f$ Kraftberechnungen
 *
 * \subsection classes_subsec Weitere modifizierte Klassen
 *
 * - World
 * - Observer
 *
 * \subsection manual_sec Anleitung
 *
 * \code
 * make && ./blatt1
 * \endcode
 *
 */
#ifndef _DEFINES_HPP
#define _DEFINES_HPP

// define the dimension of the particles
#define DIM 3
// No neighbour for parallel
#define NO_NEIGHBOUR -1
// reals in double precision
typedef double real;
// squre define
#define sqr(_x) ((_x)*(_x))

// small macro, to expand index calculation for different dimensions
#if DIM == 2
#define J(jCell,nCells) ((jCell)[1] + (nCells)[1]*(jCell)[0]) 
#elif DIM == 3
#define J(jCell,nCells) ((jCell)[2] + (nCells)[2]*((jCell[1] + (nCells)[1]*(jCell)[0])))
#endif

// another macro, to get index of ip in our process grid (by means of rank(r) and the #processes(np)
#if DIM == 2
#define Jinv(r, np, ip) \
    ((ip)[0]=(r)%(np)[0], (ip)[1]=(r)/(np)[0])
#elif DIM == 3
#define Jinv(r, np, ip) \
    ((ip)[0] = (r)%(np)[0], \
      (ip)[1] = (r)/(np)[0] % (np)[1], \
      (ip)[2] = (r)/(np)[0] / (np)[1] )
#endif


#if DIM == 3
#define Iterate(ic,minnc,maxnc) \
for( (ic)[0]=(minnc)[0]; (ic)[0]<(maxnc)[0]; (ic)[0]++ ) \
for( (ic)[1]=(minnc)[1]; (ic)[1]<(maxnc)[1]; (ic)[1]++ ) \
for( (ic)[2]=(minnc)[2]; (ic)[2]<(maxnc)[2]; (ic)[2]++ )
#elif DIM == 2
#define Iterate(ic,minnc,maxnc) \
for( (ic)[0]=(minnc)[0]; (ic)[0]<(maxnc)[0]; (ic)[0]++ ) \
for( (ic)[1]=(minnc)[1]; (ic)[1]<(maxnc)[1]; (ic)[1]++ )
#endif

#endif // _DEFINES_HPP
// vim:set et sts=4 ts=4 sw=4 ai ci cin cino=g0,t0:
