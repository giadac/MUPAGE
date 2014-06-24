/*#############################################################
 *#                                                           #
 *#    Author: G.Carminati                                    #
 *#                                                           #
 *#############################################################
 */

#include "Muons.hh"
#include "Parameters.hh"

Muons::Muons(double zenith, double azimuth, unsigned int multiplicity, 
	     double depth, TRandom3& trand)
  : zenith_(zenith), azimuth_(azimuth), multiplicity_(multiplicity), 
    depth_(depth), trand_(trand)
{}


//=============================================================
/*    First Release: Jan 2007                                 #          
 *#############################################################
 *
 *   This function sets parameters
 */
void Muons::Set(const double& zenithSG, const double& azimuth, 
		const unsigned int& multiplicity, const double& depth)
 {
   zenith_       = zenithSG;
   azimuth_      = azimuth;
   multiplicity_ = multiplicity;
   depth_        = depth;
 }   


//=============================================================
/*    First Release: Jan 2007                                 #          
 *#############################################################
 *
 *   This function sets the parameters of file parameters.inp
 */
void Muons::GetParameters(const double& depthmax, const double& canzmin, 
			  const double& canzmax, const double& canradius, 
			  const double& lsmin, const double& lsmax, 
			  const double& energymin, const double& energymax)
 {
   depthmax_  = depthmax;
   canzmin_   = canzmin;
   canzmax_   = canzmax;
   canradius_ = canradius;
   lsmin_     = lsmin;
   lsmax_     = lsmax;
   energymin_ = energymin;
   energymax_ = energymax;
 }


//=============================================================
/*    First Release: May 2005                                 #          
 *#############################################################
 *
 *   This function computes the flux for multiple muons
 * 
 *   INPUT:  muon vertical depth    (in km.w.e.)
 *           zenithal angle         (in rad)
 *           muon multiplicity
 *
 *   OUTPUT: flux
 */
double Muons::Flux()
{
  double ct = cos(zenith_);
  /*compute parameter kappa*/
  double k0 = pow(depth_,K0b)*K0a;
  double k1 = (K1a*depth_ + K1b );
  double k  = k0*ct*exp(k1/ct);
  /*compute parameter ni*/
  double ni0 = ni0a*depth_*depth_ + ni0b*depth_ + ni0c;
  double ni1 = ni1a*exp(ni1b*depth_);
  double ni  = ni0*exp(ni1/ct);
  /*compute the flux*/
  double flux = k / pow(multiplicity_,ni);
  return flux;
}


//=============================================================
/*    First Release: Jun 2007                                 #          
 *#############################################################
 *
 *   This function computes the projected area
 * 
 *   INPUT:  CAN radius       (in m)
 *           zenithal angle   (in rad)
 *
 *   OUTPUT: projectedarea
 */
double Muons::ProjectedArea()
{
  double area = pi*canradius_*canradius_;
  double projectedarea = area*cos(zenith_);
  return projectedarea;
}


//=============================================================
/*    First Release: Jun 2007                                 #          
 *#############################################################
 *
 *   This function computes the solid angle
 * 
 *   INPUT:  minimum angle     (in rad)
 *           maximum angle     (in rad)
 *
 *   OUTPUT: solidangle
 */
double Muons::SolidAngle(const double& alpha1, const double& alpha2)
{
  double cosines = cos(alpha1) - cos(alpha2);
  cosines = fabs(cosines); 
  double solidangle = cosines * twopi;
  return solidangle;
}


//=============================================================
/*    First Release: Jun 2007                                 #          
 *#############################################################
 *
 *   This function computes the product of the flux, the 
 *   projected area and the solid angle
 * 
 *   INPUT:  CAN radius       (in m)
 *           zenithal angle   (in rad)
 *           minimum angle     (in rad)
 *           maximum angle     (in rad)
 *
 *   OUTPUT: g
 */
double Muons::GAMMA(const double& alpha1, const double& alpha2)
{
  double g = Flux()*ProjectedArea()*SolidAngle(alpha1,alpha2);
  return g;
}


//=============================================================
/*    First Adaptation: Sep 2005                              #
 *    from ANTARES Software (canToDraw.cc)                    #
 *#############################################################
 *
 *     This function generates a random point on a 
 *     cylindric surface (not on the lower disk)
 *  
 *     INPUT:  zenithal angle                     (in rad)
 *             azimuthal angle                    (in rad)
 *
 *     OUTPUT: posBundleOnCAN = (x,y,z) on the cylindric 
 *                              surface             (in m)
 */
Vec3D Muons::axisBundleOnCAN()
{
  double Xr = 0.;
  double Yr = 0.;
  double sign = 0.;
  double Yoff = 0.;
  double phi_in = 0.;

  double ct = cos(zenith_);
  double st = sin(zenith_);
  double     a = canradius_ * ct;                
  double     b = canradius_;                     
  double Hproj = (canzmax_ - canzmin_)*0.5 * st;

  double BundleOnCAN_posx = 0.;
  double BundleOnCAN_posy = 0.;
  double BundleOnCAN_posz = 0.;

  for (bool loop = true; loop;)
    {
      Xr = trand_.Rndm() *2*canradius_ - canradius_;
      Yr = trand_.Rndm() *2*(Hproj+a) - (Hproj+a);
      if(fabs(Yr) > Hproj)                     
	{
	  sign = fabs(Yr)/Yr;
	  Yoff = Yr - sign*Hproj;
	  if( Xr*Xr/(b*b) + Yoff*Yoff/(a*a) > 1. ) continue;
	} 
      loop = false;
    }
  
  /* coordinates on the upper disk */
  if(( Yr > (Hproj-a) )&&( Xr*Xr/(b*b)+(Yr-Hproj)*(Yr-Hproj)/(a*a) <= 1. ))
    {
      /* we change sign of x coordinate to use frames of same helicity */
      BundleOnCAN_posx = - Xr;
      BundleOnCAN_posy = (Yr-Hproj)/ct;
      BundleOnCAN_posz = canzmax_;
    }
  else  /* coordinates on the lateral surface */
    {
      /* we change sign of x coordinate to use frames of same helicity */
      phi_in = azimuth_ + 3*halfpi - acos(-Xr/canradius_);
      BundleOnCAN_posx = canradius_ * cos(phi_in);
      BundleOnCAN_posy = canradius_ * sin(phi_in);
      BundleOnCAN_posz = ( Yr + a*sqrt(1.-Xr*Xr/(b*b)) )/st;
      BundleOnCAN_posz += centerZ();
    }

  Vec3D posBundleOnCAN(BundleOnCAN_posx, BundleOnCAN_posy, BundleOnCAN_posz);

  return posBundleOnCAN;
}


//=============================================================
/*    First Release: Dec 2006                                 #          
 *#############################################################
 *
 *   This function inizializes the vertical depth of muon with 
 *   respect to the sea level (in km.w.e.)
 */
double Muons::Depth(const double& zBundle) 
{
  depth_ = depthmax_ + canzmin_*MeterToKm() - zBundle*MeterToKm();
  depth_ *= kmTOkmwe();
  return depth_;
}


//=============================================================
/*    First Release: Jun 2005                                 #
 *#############################################################
 *
 *   This function computes the depth for energy distribution 
 *
 *   INPUT:  muon vertical depth   (in km.w.e.)
 *           zenithal angle        (in rad)
 *
 *   OUTPUT: betachi
 */
double Muons::BetaChi()
{
  double chi = depth_/cos(zenith_);
  double betachi = BETA*chi;
  return betachi;
}


//=============================================================
/*    First Release: Jun 2005                                 #
 *#############################################################
 *
 *   This function computes the gamma parameter for energy 
 *   distribution of single muons 
 *
 *   INPUT:  muon vertical depth  (in km.w.e.)
 *
 *   OUTPUT: gamma_singlemu
 */
double Muons::Gamma_singlemu()
{
  double gamma_singlemu = g0*log(depth_) + g1;
  return gamma_singlemu;
}


//=============================================================
/*    First Release: Jun 2005                                 #
 *#############################################################
 *
 *   This function computes the epsilon parameter for energy 
 *   distribution of single muons
 * 
 *   INPUT:  muon vertical depth  (in km.w.e.)
 *           zenithal angle       (in rad)
 *
 *   OUTPUT: epsilon_singlemu
 */
double Muons::Epsilon_singlemu()
{
  double e0 = e0a * exp(e0b*depth_);
  double e1 = e1a*depth_ + e1b;
  double epsilon_singlemu = e0/cos(zenith_) + e1;
  return epsilon_singlemu;
}


//=============================================================
/*    First Release: Jun 2005                                 #
 *#############################################################
 *
 *   This function computes the maximum energy value 
 *   of a single muon
 *
 *   OUTPUT: esinglemax
 */
double Muons::E_singlemu()
{
  double funzi = 1 - exp(-BetaChi());
  funzi *= Epsilon_singlemu();
  double esinglemax = funzi / (Gamma_singlemu() - 1);
  return esinglemax;
}


//=============================================================
/*    First Release: Jun 2005                                 #
 *#############################################################
 *
 *   This function computes the energy distribution of a
 *   single muon
 *
 *   INPUT:  E_mu = energy value        (in TeV)
 *
 *   OUTPUT: dnde_singlemu
 */
double Muons::dNdE_singlemu(const double& E_mu)
{
  /*compute the normalization factor G*/
  double gamma = Gamma_singlemu();
  double epsilon = Epsilon_singlemu();
  double bx = BetaChi();
  double funzi = 1 - exp(-bx);
  double gammaLess1 = gamma - 1;
  double func_g = pow(funzi, gammaLess1);
  double g = 2.3*gammaLess1*pow(epsilon, gammaLess1);
  g *= exp(gammaLess1*bx)*func_g;
  /*compute the maximum value of distribution energy*/
  double funzi2 = funzi*epsilon + E_mu;
  double func = pow(funzi2, -gamma);
  double espo = exp(bx*(1-gamma)) * E_mu;
  double dnde_singlemu = g*espo*func;
  return dnde_singlemu;
}


//=============================================================
/*    First Release: Jul 2005                                 #
 *#############################################################
 *
 *   This function computes the energy of a single muon
 *
 *   OUTPUT: energy                  (in TeV)
 */
double Muons::HoR_Energy_singlemu()
{
  const double logEmin = log10(energymin_);
  const double logEmax = log10(energymax_);
  double energy;

  /* Hit or Miss method to accept (or reject) the energy value computed */
  bool loop = true;
  while(loop)
    {
      double loge = trand_.Rndm() * (logEmax - logEmin) + logEmin;
      energy = pow(10, loge);
      double f_e = dNdE_singlemu(energy);
      double u_e = trand_.Rndm() * dNdE_singlemu( E_singlemu() );
      if (u_e <= f_e)
	loop = false;
    }
  return energy;
}


//=============================================================
/*    First Release: Jul 2005                                 #
 *#############################################################
 *
 *   This function computes the average radius value for  
 *   radial distribution  
 *
 *   INPUT:  multiplicity
 *           zenithal angle       (in rad)
 *           muon vertical depth  (in km.w.e.)
 *   OUTPUT: average_r
 */
double Muons::AverageRadius()
{
  int mult = ( multiplicity_ < 4 ) ? multiplicity_ : 4;
  double rho0 = RHO0a * mult + RHO0b;
  double espo = exp((zenith_ - THETA0) * Fpiccolo) + 1;
  double F = 1 / espo;
  double average_r = rho0 * pow(depth_,RHO1) * F;
  return average_r;
}


//=============================================================
/*    First Release: Jul 2005                                 #
 *#############################################################
 *
 *   This function computes the alpha parameter for  
 *   radial distribution  
 *
 *   INPUT:  muon vertical depth  (in km.w.e.)
 *           multiplicity
 *
 *   OUTPUT: alpha
 */
double Muons::Alpha()
{
  int mult = ( multiplicity_ < 4 ) ? multiplicity_ : 4;
  double alpha0 = ALPHA0a*mult + ALPHA0b;
  double alpha1 = ALPHA1a*mult + ALPHA1b;
  double alpha  = alpha0 * exp(alpha1*depth_);
  return alpha;
}


//=============================================================
/*    First Release: Jul 2005                                #
 *#############################################################
 *
 *   This function computes the maximum radial distribution 
 *   value of a muon in a bundle
 *
 *   OUTPUT: rpeak
 */
double Muons::Rpeak()
{
  double alpha = Alpha();
  double R0 = (alpha - 3) * AverageRadius() * 0.5;
  double rpeak = R0 / (alpha - 1);
  return rpeak;
}


//=============================================================
/*    First Release: Jul 2005                                 #
 *#############################################################
 *
 *   This function computes the radial distribution of a muon 
 *   in a bundle
 *
 *   INPUT:  r = lateral spread                (in meters)
 *
 *   OUTPUT: dndr
 */
double Muons::dNdR(const double& r)
{
  double alpha = Alpha();
  double R0 = (alpha - 3) * AverageRadius() * 0.5;
  double AlphaLessTwo = alpha - 2;
  double C = (alpha - 1) * AlphaLessTwo * pow(R0, AlphaLessTwo);
  double denom = pow((r + R0), alpha);
  double dndr = C * r / denom;
  return dndr;
}


//=============================================================
/*    First Release: Jul 2005                                 #
 *#############################################################
 *
 *   This function generates the lateral spread (in function
 *   of the radial distribution) of a muon in bundle
 *
 *   OUTPUT: radius
 */
double Muons::HoR_Radius()
{
  double radius;

  /* Hit or Miss method to accept (or reject) the radial value computed */
  bool loop = true;
  while(loop)
    {
      radius = trand_.Rndm() * (lsmax_ - lsmin_) + lsmin_;
      double f_R = dNdR(radius);
      double u_R = trand_.Rndm() * dNdR( Rpeak() );
      if (u_R <= f_R)
	loop = false;
    }
  return radius;
}


//=============================================================
/*    First Release: Jun 2005                                 #
 *#############################################################
 *
 *   This function computes the gamma parameter for energy 
 *   distribution of multiple muons 
 *
 *   INPUT:  r = lateral spread    (in meter)
 *           muon vertical depth  (in km.w.e.)
 *           multiplicity
 *
 *   OUTPUT: gamma_multimu
 */
double Muons::Gamma_multimu(const double& r)
{
  int mult = ( multiplicity_ < 4 ) ? multiplicity_ : 4;
  double a  = A0 * depth_ + A1;
  double b0 = B0a * mult + B0b;
  double b1 = B1a * mult + B1b;
  double b  = b0 * depth_ + b1;
  double q  = Q0 * depth_ + Q1;
  double func_gamma = 1 - (exp(q*r)*.5);
  double gamma_multimu = a*r + b*func_gamma;
  return gamma_multimu;
}


//=============================================================
/*    First Release: Jun 2005                                 #
 *#############################################################
 *
 *   This function computes the epsilon parameter for energy 
 *   distribution of single muons 
 *
 *   INPUT:  r = lateral spread      (in meter)
 *           muon vertical depth  (in km.w.e.)
 *           multiplicity
 *
 *   OUTPUT: epsilon_multimu
 */
double Muons::Epsilon_multimu(const double& r)
{
  double c0 = C0a*depth_ + C0b;
  double c  = c0 * exp(C1*r);
  double d0 = D0a*depth_ + D0b;
  double d1 = D1a*depth_ + D1b;
  double d  = d0 * pow(r, d1);
  double epsilon_multimu = c*zenith_ + d;
  return epsilon_multimu;
}


//=============================================================
/*    First Release: Jun 2005                                 #
 *#############################################################
 *
 *   This function computes the maximum energy value of a
 *   multiple muon
 *
 *   INPUT:  r = lateral spread    (in meter)
 *
 *   OUTPUT: emultimax
 */
double Muons::E_multimu(const double& r)
{
  double funzi = 1 - exp(-BetaChi());
  funzi *= Epsilon_multimu(r);
  double emultimax = funzi / (Gamma_multimu(r) - 1);
  return emultimax;
}


//=============================================================
/*    First Release: Jun 2005                                 #
 *#############################################################
 *
 *   This function computes the energy distribution of a
 *   multi muon in bundle
 *
 *   INPUT:  r    = lateral spread   (in meter)
 *           E_mu = energy value     (in TeV)
 *
 *   OUTPUT: dnde_multimu
 */
double Muons::dNdE_multimu(const double& r, const double& E_mu)
{
  /*compute the normalization factor G*/
  double gamma = Gamma_multimu(r);
  double epsilon = Epsilon_multimu(r);
  double bx = BetaChi();
  double GammaLessOne = gamma - 1;
  double funzi = 1 - exp(-bx);
  double func_g = pow(funzi, GammaLessOne);
  double g = 2.3*GammaLessOne*pow(epsilon, GammaLessOne);
  g *= exp(GammaLessOne*bx)*func_g;
  /*compute the maximum value of distribution energy*/
  double funzi2 = funzi * epsilon + E_mu;
  double func = pow(funzi2, -gamma);
  double espo = exp(bx*(1-gamma)) * E_mu;
  double dnde_multimu = g * espo * func;
  return dnde_multimu;
}


//=============================================================
/*    First Release: Jul 2005                                 #
 *#############################################################
 *
 *   This function computes the energy of a multiple muon 
 *   in a bundle
 *
 *   INPUT:  r = lateral spread                  (in meter)
 *
 *   OUTPUT: energy (in TeV)
 */
double Muons::HoR_Energy_multimu(const double& r)
{
  const double logEmin = log10(energymin_);
  const double logEmax = log10(energymax_);
  double energy;

  /* Hit or Miss method to accept (or reject) the energy value computed */
  bool loop = true;
  while(loop)
    {
      double loge = trand_.Rndm() * (logEmax - logEmin) + logEmin;
      energy = pow(10, loge);
      double f_e = dNdE_multimu(r, energy);
      double u_e = trand_.Rndm() * dNdE_multimu(r, E_multimu(r) );
      if (u_e <= f_e)
	loop = false;
    }
  return energy;
}


//=============================================================
/*    First Release: Mar 2006                                 #
 *#############################################################
 *
 *   This function computes the position of each muon in the 
 *   bundle on the laboratory frame
 *
 *   INPUT:  posBundleOnCAN = coordinates (x,y,z) of the 
 *                               shower axis impact point on 
 *                               the can surface        (in m)
 *           zenithal angle                           (in rad)
 *           azimuthal angle                          (in rad)
 *           r = lateral spread                         (in m)
 *
 *   OUTPUT: posMuOnLAB = coordinates (x,y,z) of each muon
 *                        on the laboratory frame       (in m)
 */
Vec3D Muons::multimuOnLAB(Vec3D posBundleOnCAN, const double& r)
{
  /* Xpi e Ypi are the coordinates of the intercept of the
     track with the orthogonal plane to bundle axis */
  double beta = trand_.Rndm() * twopi;    /* in rad */
  double Xpi  = r*cos(beta);              /* in m   */
  double Ypi  = r*sin(beta);              /* in m   */

  /* introduce the Euler angles (Phi, Theta, Psi) */
  double Phi   = - halfpi;
  double Theta = pi - zenith_;
  double Psi   = azimuth_ + halfpi;

  double cf = cos(Phi);
  double sf = sin(Phi);
  double ct = cos(Theta);
  double st = sin(Theta);
  double cp = cos(Psi);
  double sp = sin(Psi);

  /* introduce the elements of the rotational matrix A, which is the
     composition of the three rotations with Euler angles in the 
     so-called 'X-convention' */
  double Axx =   cf*cp - ct*sf*sp;
  double Axy = - sf*cp - ct*cf*sp;
//   double Axz =   st*sp;         <------- unused variable
  
  double Ayx =   cf*sp + ct*sf*cp;
  double Ayy = - sf*sp + ct*cf*cp;
//   double Ayz = - st*cp;         <------- unused variable
  
  double Azx = st*sf;
  double Azy = st*cf;
//   double Azz = ct;              <------- unused variable

  double posx_MuOnLAB = Axx*Xpi + Axy*Ypi + posBundleOnCAN.x;  /* in m */
  double posy_MuOnLAB = Ayx*Xpi + Ayy*Ypi + posBundleOnCAN.y;  /* in m */
  double posz_MuOnLAB = Azx*Xpi + Azy*Ypi + posBundleOnCAN.z;  /* in m */

  Vec3D posMuOnLAB(posx_MuOnLAB, posy_MuOnLAB, posz_MuOnLAB);

  return posMuOnLAB;
}


//=============================================================
/*    First Adaptation: Sep 2005                              #
 *    from ANTARES Software                                   #
 *#############################################################
 *
 *   This function computes the position of the muon (single 
 *   or one in the bundle) on the cylindric surface 
 *   (not on the lower disk)
 *
 *   INPUT:  posMuBefore = coordinates (x,y,z) of the muon on 
 *                         the laboratory frame          (in m)
 *           uBundle     = direction cosines of bundle axis
 *
 *   OUTPUT: posMuAfter = coordinates (x,y,z) of muon impact 
 *                        point on the enlarged can surface
 *                                                       (in m)
 *
 * @return error code
 *  <li> 0   - successfull operation
 *  <li> -99 - the muon does not intercept the can
 * @
 */
int Muons::muOnCAN(Vec3D posMuBefore, Vec3D uBundle, Vec3D& posMuAfter)
{
  /* <----------- position on the upper disk -----------> */
  /* intersect the line (x_i,y_i,z_i) = (posMuBefore) + l*uBundle with
     the plane z = canzmax */
  double l = (canzmax_ - posMuBefore.z)/uBundle.z;
  double xx = posMuBefore.x + l*uBundle.x;   /* in m */
  double yy = posMuBefore.y + l*uBundle.y;   /* in m */
  if (xx*xx + yy*yy <= canradius_*canradius_)
    {
      posMuAfter.x = xx;
      posMuAfter.y = yy;
      posMuAfter.z = canzmax_;
    }
  /* <----------- position on barrel surface -----------> */
  else
    {
      posMuBefore.z -= centerZ();
      /* intersect the line (x_i,y_i,z_i) = (posMuBefore) + k*uBundle with
	 the cylinder (of equation xx**2 + yy**2 = R**2)
	 => solving the second degree equation a*k**2 + 2*b*k + c = 0; 
	 the choice between the 2 solutions is made using the 
	 fact that the angle between u and the solution must be < PI/2 */
      double a =     uBundle.x*uBundle.x     +     uBundle.y*uBundle.y;
      double b =     uBundle.x*posMuBefore.x +     uBundle.y*posMuBefore.y;
      double c = posMuBefore.x*posMuBefore.x + posMuBefore.y*posMuBefore.y 
	         - canradius_*canradius_;
      double delta = b*b - a*c;
      if (delta < 0)
	return -99;
      double sqrtdelta = sqrt(delta);
  
      double  k1 = ( -b - sqrtdelta )/a;
      double ax1 = posMuBefore.x + k1*uBundle.x;
      double ay1 = posMuBefore.y + k1*uBundle.y;
      double az1 = posMuBefore.z + k1*uBundle.z;
      
      double  k2 = ( -b + sqrtdelta )/a;
      double ax2 = posMuBefore.x + k2*uBundle.x;
      double ay2 = posMuBefore.y + k2*uBundle.y;
      double az2 = posMuBefore.z + k2*uBundle.z;
      
      if(az1 > az2)
	{
	  posMuAfter.x = ax1;
	  posMuAfter.y = ay1;
	  posMuAfter.z = az1;
	}
      else if(az2 > az1)
	{
	  posMuAfter.x = ax2;
	  posMuAfter.y = ay2;
	  posMuAfter.z = az2;
	}
      else
	{ return -99; }
      posMuAfter.z += centerZ();
    }
  if (posMuAfter.z > canzmax_ || posMuAfter.z < canzmin_)
    return -99;
  return 0;
}

