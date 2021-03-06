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
// reals in double precision
typedef double real;
// squre define
#define sqr(_x) ((_x)*(_x))

// small macro, to expand index calculation for different dimensions
#if DIM == 2
#define J(jCell,cell_N) ((jCell)[1] + (cell_N)[1]*(jCell)[0]) 
#elif DIM == 3
#define J(jCell,cell_N) ((jCell)[2] + (cell_N)[2]*((jCell[1] + (cell_N)[1]*(jCell)[0])))
#endif

#endif // _DEFINES_HPP
// vim:set et sts=4 ts=4 sw=4 ai ci cin cino=g0,t0:
