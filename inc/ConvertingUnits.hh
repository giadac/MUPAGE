//#############################################################
//#                                                           #
//#    Author: G.Carminati                                    #
//#    First Release: Jan 2006                                #
//#                                                           #
//#############################################################
//

#ifndef CONVERTING_UNITS_H
#define CONVERTING_UNITS_H

#include <cmath>

#define pi M_PI
#define halfpi M_PI_2
#define twopi 2*M_PI
#define c_light 2.99792458e+8

/* Length [L] */
inline double MeterToKm()   { return 1.e-3; }
inline double KmToMeter()   { return 1.e+3; }
inline double Km2ToMeter2() { return KmToMeter()*KmToMeter(); }

/* Angle */
inline double RadToDeg() { return 180./pi; }
inline double DegToRad() { return pi/180.; }

/* Time [T] */
inline double sTOdays() { return 1/86400.; }

/* Energy [E] */
inline double GeVtoTeV() { return 1.e-3; }
inline double TeVtoGeV() { return 1.e+3; }

#endif /* CONVERTING_UNITS_H */
