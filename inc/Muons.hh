/*#############################################################
 *#                                                           #
 *#    Author: G.Carminati                                    #
 *#    First Release: Dec 2006                                #
 *#                                                           #
 *#############################################################
 */

#ifndef MUONS_h
#define MUONS_h

#include "TRandom3.h"

#include "ConvertingUnits.hh"
#include "Geometry.hh"

class Muons
{
private: 
  double       zenith_;
  double       azimuth_;
  unsigned int multiplicity_;
  double       depth_;
  TRandom3&    trand_;

  double depthmax_;
  double canzmin_;
  double canzmax_;
  double canradius_;
  double lsmin_;
  double lsmax_;
  double energymin_;
  double energymax_;

public:
  Muons(double zenith, double azimuth, unsigned int multiplicity, double depth,
	TRandom3& trand);
  virtual ~Muons(){};

  void   Set(const double& zenith, const double& azimuth, 
	     const unsigned int& multiplicity, const double& depth);
  void   GetParameters(const double& depthmax, const double& canzmin, 
		       const double& canzmax, const double& canradius, 
		       const double& lsmin, const double& lsmax, 
		       const double& energymin, const double& energymax);
  double Flux();
  double ProjectedArea();
  double SolidAngle(const double& alpha1, const double& alpha2);
  double GAMMA(const double& alpha1, const double& alpha2);
  Vec3D  axisBundleOnCAN();
  double Depth(const double& zBundle); 
  double BetaChi();
  double Gamma_singlemu();
  double Epsilon_singlemu();
  double Gamma_multimu(const double& r);
  double Epsilon_multimu(const double& r);
  double E_singlemu();
  double dNdE_singlemu(const double& E_mu);
  double HoR_Energy_singlemu();
  double E_multimu(const double& r);
  double dNdE_multimu(const double& r, const double& E_mu);
  double HoR_Energy_multimu(const double& r);
  double AverageRadius();
  double Alpha();
  double Rpeak();
  double dNdR(const double& r);
  double HoR_Radius();
  Vec3D  multimuOnLAB(Vec3D posBundleOnCAN, const double& r);
  int    muOnCAN(Vec3D posMuBefore, Vec3D uBundle, Vec3D& posMuAfter);
};

#endif /*MUONS_h*/
