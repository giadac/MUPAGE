//#############################################################
//#                                                           #
//#    Author: G.Carminati                                    #
//#    First Release: Feb 2006                                #
//#                                                           #
//#############################################################
//

#ifndef PARAMETERS_H
#define PARAMETERS_H

#include<string>

#include "ConvertingUnits.hh"

typedef double ENERGY;

static const ENERGY Emin = 0.001f;
static const ENERGY Emax = 500.0f;

typedef double ANTADouble;
typedef unsigned short int ANTAUsInt;
typedef unsigned long int ANTAUlInt;
typedef std::string ANTAString;

typedef double FLUX;
typedef double ALLmu;
typedef double SINGLEmu;
typedef double MULTImu;
typedef double LateralSpread;
 
extern ANTADouble DEPTHmax;
extern ANTADouble Zmin;
extern ANTADouble Zmax;
extern ANTADouble CANr;
extern ANTADouble EnlargedCANr;
extern ANTADouble CANh;
extern ANTADouble THETAmin;
extern ANTADouble THETAmax;
extern ANTADouble Rmin;
extern ANTADouble Rmax;
extern ANTAUsInt  MULTmin;
extern ANTAUsInt  MULTmax;
extern ANTAUsInt  GEANTid;
extern ANTADouble density;
extern ANTADouble AbsLength;
extern ANTADouble NAbsLength;
extern ANTADouble Ethreshold;

extern ANTAUlInt  seed;
extern ANTAUlInt  events;
extern ANTAString numrun;
extern ANTAString OutputFileName;
extern ANTAString LivetimeFileName;
extern ANTAString ParameterFileName;

//WARNING: Do not change the units!

/* to convert km in km.w.e. (1 km.w.e. = 1.e+6 kg m**(-3)) */
inline double kmTOkmwe()   { return 1.*density; }

inline double CANheight()  { return Zmax - Zmin; }                  /* in m  */
inline double centerZ()    { return (Zmax + Zmin)*0.5; }            /* in m  */
inline double CANradius()  { return CANr*MeterToKm(); }             /* in km */
inline double DEPTHmin()   { return DEPTHmax - CANheight()*MeterToKm(); }
                                                                    /* in km */
inline double Qmin()       { return THETAmin*DegToRad(); }         /* in rad */
inline double Qmax()       { return THETAmax*DegToRad(); }         /* in rad */
inline double weDEPTHmin() { return DEPTHmin()*kmTOkmwe(); }   /* in km.w.e. */
inline double weDEPTHmax() { return DEPTHmax*kmTOkmwe(); }     /* in km.w.e. */
inline double max(const double& a, const double& b) { return a > b ? a : b; }
inline double min(const double& a, const double& b) { return a < b ? a : b; }


/* Parameters from APP 25 (2006) 1-13 */
static const FLUX  K0a =  0.0072f;
static const FLUX  K0b = -1.927f;
static const FLUX  K1a = -0.581f;
static const FLUX  K1b =  0.034f;
static const FLUX ni0a = -0.0771f;
static const FLUX ni0b =  0.524f;
static const FLUX ni0c =  2.068f;
static const FLUX ni1a =  0.03f;
static const FLUX ni1b =  0.47f; 

static const ALLmu BETA = 0.42f;

static const SINGLEmu  g0 = -0.232f;
static const SINGLEmu  g1 =  3.961f;
static const SINGLEmu e0a =  0.0304f;
static const SINGLEmu e0b =  0.359f;
static const SINGLEmu e1a = -0.0077f;
static const SINGLEmu e1b =  0.659f;

static const MULTImu  A0 =  0.0033f;
static const MULTImu  A1 =  0.0079f;
static const MULTImu B0a =  0.0407f;
static const MULTImu B0b =  0.0283f;
static const MULTImu B1a = -0.312f;
static const MULTImu B1b =  6.124f;
static const MULTImu  Q0 =  0.0543f;
static const MULTImu  Q1 = -0.365f;
static const MULTImu C0a = -0.069f;
static const MULTImu C0b =  0.488f;
static const MULTImu  C1 = -0.117f;
static const MULTImu D0a = -0.398f;
static const MULTImu D0b =  3.955f;
static const MULTImu D1a =  0.012f;
static const MULTImu D1b = -0.35f;

static const LateralSpread    RHO0a = - 1.786f;
static const LateralSpread    RHO0b =  28.26f;
static const LateralSpread     RHO1 = - 1.06f;
static const LateralSpread   THETA0 =   1.3f;
static const LateralSpread Fpiccolo =  10.4f;
static const LateralSpread  ALPHA0a = - 0.448f;
static const LateralSpread  ALPHA0b =   4.969f;
static const LateralSpread  ALPHA1a =   0.0194f;
static const LateralSpread  ALPHA1b =   0.276f;

#endif /*PARAMETERS_H*/
