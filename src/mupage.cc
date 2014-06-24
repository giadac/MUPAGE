/*
 *MM        MM UU        UU PPPPPPPPPPP   AAAAAAAAAA   GGGGGGGGGG  EEEEEEEEEEEE
 *MMM      MMM UU        UU PPPPPPPPPPPP AAAAAAAAAAAA GGGGGGGGGGGG EEEEEEEEEEEE
 *MMMM    MMMM UU        UU PP        PP AA        AA GG        GG EE
 *MM MM  MM MM UU        UU PP        PP AA        AA GG           EE
 *MM  MMMM  MM UU        UU PP        PP AA        AA GG           EE
 *MM   MM   MM UU        UU PPPPPPPPPPPP AAAAAAAAAAAA GG           EEEEEEEE
 *MM        MM UU        UU PPPPPPPPPPP  AAAAAAAAAAAA GG     GGGGG EEEEEEEE
 *MM        MM UU        UU PP           AA        AA GG     GGGGG EE
 *MM        MM UU        UU PP           AA        AA GG        GG EE
 *MM        MM UU        UU PP           AA        AA GG        GG EE
 *MM        MM UUUUUUUUUUUU PP           AA        AA GGGGGGGGGGGG EEEEEEEEEEEE
 *MM        MM  UUUUUUUUUU  PP           AA        AA  GGGGGGGGGG  EEEEEEEEEEEE
 *
 *@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
 * 
 * Author:  G.Carminati
 * version 1:  3rd December 2007
 * version 2:  3rd June 2009
 * 
 */

#include<iostream>
#include<fstream>
#include<cmath>
#include<iomanip>
#include<string>
#include<sstream>

#include "Parameters.hh"
#include "Muons.hh"
#include "Decode.hh"

#include "TRandom3.h"

int main(int argc, char** argv)
{
  /* information for execution */
  Decode dec; 
  int nerr = dec.CommandLine(argc, argv);
  if(nerr < 0)
    {
      std::cout << dec.ErrorCommandLine();
      exit(3);
    }
  int nerror = 0;
  nerror = dec.DecodeParameters();
  if(nerror > 0)
    {
      std::cout << dec.ErrorDecoder(nerror);
      exit(2);
    }
  else if(nerror < 0)
    {
      std::cout << dec.ErrorOpenFiles(nerror);
      exit(1);
    }

  /* opening output files */
  ofstream fout, livetime_fout;
  fout.open(OutputFileName.c_str(), std::ios::out);
  livetime_fout.open(LivetimeFileName.c_str(), std::ios::out);
  if (fout.fail() )
    {
      nerror = -98; 
      std::cout << dec.ErrorOpenFiles(nerror);
      exit(1);
    }
  if (livetime_fout.fail() )
    {
      nerror = -97; 
      std::cout << dec.ErrorOpenFiles(nerror);
      exit(1);
    }

  /* seeds definition */
  if(seed == 0)
    seed = time(NULL);
  TRandom3 genphi(seed);
  TRandom3 genm(seed++);
  TRandom3 gentheta(seed+2);
  TRandom3 genu(seed+3);
  TRandom3 gen(seed+4);

  /* create parameters for the CAN */
  Vec3D can(Zmin, Zmax, CANr);
  CANh = CANheight();

  /* the CAN have to be enlarged to avoid losing events */
  Vec3D enlargedCAN(can.x, can.y, can.z+EnlargedCANr);

  /* check on vertical depth values */
  if (weDEPTHmin() < 1.5 || weDEPTHmin() > 5.)
    {
      nerror = 6;
      std::cout << dec.ErrorDecoder(nerror);
      exit(2);
    }
  if (weDEPTHmax() < 1.5 || weDEPTHmax() > 5.)
    {
      nerror = 7;
      std::cout << dec.ErrorDecoder(nerror);
      exit(2);
    }

  unsigned int nrun = atoi(&numrun[0]);
  unsigned int position = 0;
  
  fout.setf(std::ios::showpoint|std::ios::fixed);
  fout <<"Run Number: " <<nrun<<std::endl;
  fout <<"livetime : ";
  position = fout.tellp();
  fout << "                              " << std::endl;

  /* Implementation class Muons */
  Muons muon (Qmin(), 0., MULTmin, weDEPTHmin(), gen);
  muon.GetParameters(DEPTHmax, enlargedCAN.x, enlargedCAN.y, 
		     enlargedCAN.z, Rmin, Rmax, Emin, Emax);

  /* initialization of the counter to compute the generated events 
     on the enlarged can */
  unsigned int counter = 0;

  /* initialization of counters to compute livetime */
  unsigned int thetarange = 86;
  unsigned int ncounter[thetarange];
  memset(ncounter, 0, thetarange*sizeof(int));

  /* compute the area of the enlarged CAN (in m**2) */
  double Sz  = pi*enlargedCAN.z*enlargedCAN.z;      //<------ upper disk
  double Sxy = 2*enlargedCAN.z*CANh;                //<------ barrel surface

  /* compute the maximum projected area using the fact that 
     d(Area)/d(theta) = 0 (in m**2) */
  double ALFA = atan( - enlargedCAN.z /( 2*CANh ))*.5 + halfpi ;

  /* compute the maximum number of events per unit of time */
  double halfdegree = 0.5 * DegToRad();
  double Smax = Sz*cos(ALFA) + Sxy*sin(ALFA);
  Smax = max(Sz, Smax);
  double Nmax = muon.Flux()*Smax*muon.SolidAngle(ALFA-halfdegree,
						 ALFA+halfdegree);

  /* counter for events with energy lower than threshold one */
  unsigned int Nlow = 0;

  for (unsigned int i = 0; i < events; ++i)
    {
      muon.GetParameters(DEPTHmax, enlargedCAN.x, enlargedCAN.y, 
			 enlargedCAN.z, Rmin, Rmax, Emin, Emax);

      unsigned int m = 0;
      double phi     = 0.;
      /* thetaSG means the zenith angle in the System of Generation */
      double thetaSG = 0.;
      double ctSG    = 0.;
      double stSG    = 0.;
      /* thetaDF means the zenith angle in the Direction-to-fly System */
      double thetaDF = 0.;
      double ctDF    = 0.;
      double stDF    = 0.;
      double u       = 0.;
      double wedepthBundle = 0.;
      Vec3D posBundle(0.,0.,0.);

      /* generate random the azimuth (phi) angle [in radians] */
      phi = genphi.Rndm() * twopi;
      /* WARNING: the azimuthal angle is the angle of the pointlike
       *          source, so upward-going; instead, the muon is
       *          downward-going. So, it is necessary to convert
       *          the azimuthal angle into
       *          phi_mu = 180 degrees + azimuth
       */
      phi += pi;    /* in rad */


      bool loop = true;
      while(loop)
	{
	  /* generate random the muon multiplicity */
	  m = (int)(genm.Rndm()*(MULTmax + 1 - MULTmin)) + MULTmin;

	  /* generate random the zenith (theta) angle [in radians] */
 	  thetaSG = gentheta.Rndm() * ( Qmax() - Qmin() ) + Qmin();
	  ctSG = cos(thetaSG);
	  stSG = sin(thetaSG);

	  /* implement new parameters in the class */
	  muon.Set(thetaSG, phi, m, wedepthBundle);

	  /* compute the shower axis position 
	     on the can surface (in m) */
	  posBundle = muon.axisBundleOnCAN();

	  /* convert the z-position of the impact point 
	     in the vertical depth with respect to the sea level 
	     (in km.w.e.) */
	  wedepthBundle = muon.Depth(posBundle.z);

	  /* compute the projected area of the enlarged can 
	     for the generated zenithal angle (in m**2) */
	  double S = Sz*ctSG + Sxy*stSG;

	  /* compute the number of events per unit of time 
	     for generated zenithal angle, multiplicity and vertical depth */
	  double Nproj = muon.Flux()*S*muon.SolidAngle(thetaSG+halfdegree,
						       thetaSG-halfdegree);

	  /* accept or reject event */
	  u = genu.Rndm() * Nmax;
	  if (u <= Nproj)         //<------- accept event
	    loop = false;
	}

      /* compute the generated single events N* 
	 on upper disk of the enlarged CAN at different zenithal angles */
      float fhBundle = wedepthBundle;
      float    fHmin = weDEPTHmin();
      if (fhBundle == fHmin && m == MULTmin)
	{
	  ++counter;
	  for(unsigned int j = 0; j < thetarange; ++j)
	    {
	      if(thetaSG * RadToDeg() > j && thetaSG * RadToDeg() <= j+1)
		++(ncounter[j]);
	    }
	}

      /* WARNING: Accordling with the switching of the azimuthal angle, 
       *          the zenithal angle must be switched into
       *          theta_mu = 180 degrees - zenith
       */
      thetaDF = pi - thetaSG;     /* in rad */
      ctDF = cos(thetaDF);
      stDF = sin(thetaDF);
      muon.Set(thetaSG, phi, m, wedepthBundle);

      /* unit vector (i.e. direction cosines) of Bundle Axis */
      double cf = cos(phi);
      double sf = sin(phi);
      Vec3D uBundle(stDF*cf, stDF*sf, ctDF);

      if (m == 1) /* single muons */
      	{
	  const double t = 0.;
	  /* compute the energy of single muon */
	  double e = muon.HoR_Energy_singlemu();    /* in TeV */

	  /* check if the generated energy is greater than the 
	     energy of threshold */
	  if (e >= Ethreshold)
	    {
	      fout.setf(std::ios::showpoint|std::ios::fixed);
	      fout <<std::setprecision(3)
		   <<std::setw(10)<<i+1<<std::setw(5)<<m<<std::setw(5)<<m
		   <<std::setw(10)<<posBundle.x
		   <<std::setw(10)<<posBundle.y
		   <<std::setw(10)<<posBundle.z
		   <<std::setprecision(6)
		   <<std::setw(11)<<uBundle.x<<std::setw(11)<<uBundle.y
		   <<std::setw(11)<<-uBundle.z
		   <<std::setprecision(3)
		   <<std::setw(10)<<e*TeVtoGeV()
		   <<std::setw(10)<<t<<std::setw(5)<<GEANTid<<std::endl;
	    }
	  else /* events under the energy of threshold */
	    {
	      ++Nlow;
	      --i;
	      continue;
	    }
	}

      else /* multiple muons */

      	{
	  unsigned int mc = 0;
	  double bundle_energy = 0.;
	  double xMuOnLAB[m];
	  double yMuOnLAB[m];
	  double zMuOnLAB[m];
	  double xMuOnCAN[m];
	  double yMuOnCAN[m];
	  double zMuOnCAN[m];
	  double distance[m];
	  double delaytime[m];
	  double E[m];
	  memset(xMuOnLAB, 0, m*sizeof(double));
	  memset(yMuOnLAB, 0, m*sizeof(double));
	  memset(zMuOnLAB, 0, m*sizeof(double));
	  memset(xMuOnCAN, 0, m*sizeof(double));
	  memset(yMuOnCAN, 0, m*sizeof(double));
	  memset(zMuOnCAN, 0, m*sizeof(double));
	  memset(distance, 0, m*sizeof(double));
	  memset(delaytime, 0, m*sizeof(double));
	  memset(E, 0, m*sizeof(double));

      	  for (unsigned int l = 0; l < m; ++l)
      	    {
	      /* compute the lateral spread of each muon in bundle */
      	      double R[m];
	      memset(R, 0, m*sizeof(double));
	      R[l]= muon.HoR_Radius();

	      /* compute the position of each muon in bundle on the laboratory 
		 frame (in m) */
 	      Vec3D posMuOnLAB = muon.multimuOnLAB(posBundle, R[l]);
	      xMuOnLAB[l] = posMuOnLAB.x;
	      yMuOnLAB[l] = posMuOnLAB.y;
	      zMuOnLAB[l] = posMuOnLAB.z;

 	      /* compute if the i-th muon in bundle impacts point also
		 on the can surface */ 
	      Vec3D posMuOnCAN (0.,0.,0.);
	      int iok = muon.muOnCAN(posMuOnLAB, uBundle, posMuOnCAN);
	      if (iok != 0 && l == 0) /* muon does not intercept the can*/
		--l;
	      else if (iok == 0) /* muon intercepts the can */
		{
		  xMuOnCAN[l] = posMuOnCAN.x;
		  yMuOnCAN[l] = posMuOnCAN.y;
		  zMuOnCAN[l] = posMuOnCAN.z;
		
		  /* compute the delay time (in ns) */
		  distance[l]=sqrt( (xMuOnCAN[l] - xMuOnLAB[l])*
				    (xMuOnCAN[l] - xMuOnLAB[l]) +
				    (yMuOnCAN[l] - yMuOnLAB[l])*
				    (yMuOnCAN[l] - yMuOnLAB[l]) +
				    (zMuOnCAN[l] - zMuOnLAB[l])*
				    (zMuOnCAN[l] - zMuOnLAB[l]) );
		  if (zMuOnLAB[l] < zMuOnCAN[l])
		    distance[l] = - distance[l];
		  if (l == 0)
		    delaytime[l] = 0.;
		  else
		    delaytime[l] = (distance[l] - distance[0])/c_light * 
		                   1.e+9;
		      
		  /* compute the energy of multiple muons */
		  E[l] = muon.HoR_Energy_multimu(R[l]);
		  bundle_energy += E[l];

		  ++mc;
		}
	    }

	  if (bundle_energy >= Ethreshold)
	    {
	      for (unsigned int ll = 0; ll < m; ++ll)
		{
		  if((float)E[ll] != 0)
		    {
		      fout.setf(std::ios::showpoint|std::ios::fixed);
		      fout <<std::setprecision(3)
			   <<std::setw(10)<<i+1<<std::setw(5)<<m
			   <<std::setw(5)<<ll+1<<std::setprecision(3)
			   <<std::setw(10)<<xMuOnCAN[ll]
			   <<std::setw(10)<<yMuOnCAN[ll]
			   <<std::setw(10)<<zMuOnCAN[ll]
			   <<std::setprecision(6)<<std::setw(11)<<uBundle.x
			   <<std::setw(11)<<uBundle.y
			   <<std::setw(11)<<uBundle.z
			   <<std::setprecision(3)<<std::setw(10)
			   <<E[ll]*TeVtoGeV()<<std::setw(10)<<delaytime[ll]
			   <<std::setw(5)<<GEANTid<<std::endl;
		    }
		}
	    }
	  else /* events under the energy of threshold */
	    {
	      ++Nlow;
	      --i;
	      continue;
	    }
	}   
    }

  /* compute livetime */
  livetime_fout << std::endl;
  livetime_fout <<"Generated events: " << events << std::endl;
  livetime_fout <<"Generated events with multiplicity " << MULTmin 
		<< " on the upper disk of the CAN: " << counter << std::endl;
  livetime_fout << "\nEvents with energy higher than energy of threshold: "
		<< events <<std::endl;
  livetime_fout << "Events with energy lower than energy of threshold: " 
		<< Nlow << std::endl;
  livetime_fout << std::endl;
  livetime_fout << "theta (degrees)    N*   Delta N*    N*_MC     T (sec) "
		<< "  Delta T\n";
  livetime_fout << "******************************************************"
		<< "**********\n";

  /* range of 0 to 9 degrees */
  unsigned int nsmallangle = 0;
  for(unsigned int k = 0; k < 10; ++k)
    nsmallangle += ncounter[k];
  unsigned int ensmallangle = (int)(sqrt((double)nsmallangle)+0.5);
  muon.Set(5*DegToRad(), 0., MULTmin, weDEPTHmin());
  muon.GetParameters(DEPTHmax, enlargedCAN.x, enlargedCAN.y, 
		     enlargedCAN.z, Rmin, Rmax, Emin, Emax);
  double N_MC = muon.GAMMA(0.*DegToRad(),10.*DegToRad());
  double livetime10  = (double)nsmallangle/N_MC;
  double elivetime10 = (double)ensmallangle/N_MC;
  double weight10 = 0.;
  if(livetime10 != 0)
    weight10 = 1./(elivetime10*elivetime10);
  livetime_fout.setf(std::ios::fixed|std::ios::right);
  livetime_fout	<< std::setw(11) << "[0,10]" << std::setw(10) << nsmallangle 
		<< std::setw(9) << ensmallangle 
		<< std::setprecision(1) << std::setw(10) << N_MC;
  livetime_fout.unsetf(std::ios::fixed);
  livetime_fout.setf(std::ios::floatfield); 
  livetime_fout <<std::setprecision(3)<<std::setw(11)<<livetime10 
		<<std::setprecision(1)<<std::setw(10)<<elivetime10<<std::endl;

  /* range of 10 to 69 degrees */
  unsigned int  twocounter = 0;
  unsigned int etwocounter = 0;
  double livetime[thetarange];
  double elivetime[thetarange];
  double weight[thetarange];
  memset(livetime, 0, thetarange*sizeof(double));
  memset(elivetime, 0, thetarange*sizeof(double));
  memset(weight, 0, thetarange*sizeof(double));

  for(unsigned int jj = 10; jj < 70; ++jj)
    {
      double FASA = 0.;
      if(jj % 2 != 0)
	{
	  twocounter = ncounter[jj-1] + ncounter[jj];
	  etwocounter = (int)(sqrt((double)twocounter)+0.5);
	  muon.Set( (double)(jj*DegToRad()), 0., MULTmin, weDEPTHmin());
	  FASA = muon.GAMMA( (double)(jj-1)*DegToRad(),
			     (double)(jj+1)*DegToRad());
	  livetime[jj]  = (double)twocounter/FASA;
	  elivetime[jj] = (double)etwocounter/FASA;
	  weight[jj] = 1./(elivetime[jj]*elivetime[jj]);
	  livetime_fout.unsetf(std::ios::floatfield);
	  livetime_fout.setf(std::ios::fixed|std::ios::right);
	  livetime_fout << std::setw(5) << "]" <<jj-1 << ","<< jj+1 <<"]"
			<< std::setw(10) << twocounter 
			<< std::setw(9) << etwocounter 
			<< std::setprecision(1) << std::setw(10) << FASA;
	  livetime_fout.unsetf(std::ios::fixed);
	  livetime_fout.setf(std::ios::floatfield);
	  livetime_fout <<std::setprecision(3)<<std::setw(11)<<livetime[jj] 
			<<std::setprecision(1)<<std::setw(10)<<elivetime[jj] 
			<<std::endl;
	}
    }

  /* range of 70 to 75 degrees */
  unsigned int nangle70 = 0;
  for(unsigned int kk = 70; kk < 76; ++kk)
    nangle70 += ncounter[kk];
  unsigned int enangle70 = (int)(sqrt((double)nangle70)+0.5);
  muon.Set(72.5*DegToRad(), 0., MULTmin, weDEPTHmin());
  double N_MC70 = muon.GAMMA(70.*DegToRad(),75.*DegToRad());
  double livetime70  = (double)nangle70/N_MC70;
  double elivetime70 = (double)enangle70/N_MC70;
  double weight70 = 0.;
  if(livetime70 != 0)
    weight70 = 1./(elivetime70*elivetime70);
  livetime_fout.unsetf(std::ios::floatfield);
  livetime_fout.setf(std::ios::fixed|std::ios::right);
  livetime_fout << std::setw(11) << "]70,76]" << std::setw(10) << nangle70 
		<< std::setw(9) << enangle70 
		<< std::setprecision(1) << std::setw(10) << N_MC70;
  livetime_fout.unsetf(std::ios::fixed);
  livetime_fout.setf(std::ios::floatfield);
  livetime_fout	<<std::setprecision(3)<<std::setw(11)<<livetime70
		<<std::setprecision(1)<<std::setw(10)<<elivetime70<<std::endl; 

  /* range of 76 to 85 degrees */
  unsigned int nangle76 = 0;
  for(unsigned int k76 = 76; k76 < 86; ++k76)
    nangle76 += ncounter[k76];
  unsigned int enangle76 = (int)(sqrt((double)nangle76)+0.5);
  muon.Set(80.*DegToRad(), 0., MULTmin, weDEPTHmin());
  double N_MC80 = muon.GAMMA(76.*DegToRad(),85.*DegToRad());
  double livetime80  = (double)nangle76/N_MC80;
  double elivetime80 = (double)enangle76/N_MC80;
  double weight80 = 0.;
  if(livetime80 != 0)
    weight80 = 1./(elivetime80*elivetime80);
  livetime_fout.unsetf(std::ios::floatfield);
  livetime_fout.setf(std::ios::fixed|std::ios::right);
  livetime_fout << std::setw(11) << "]76,85]" << std::setw(10) << nangle76 
		<< std::setw(9) << enangle76 
		<< std::setprecision(1) << std::setw(10) << N_MC80; 
  livetime_fout.unsetf(std::ios::fixed);
  livetime_fout.setf(std::ios::floatfield);
  livetime_fout	<<std::setprecision(3)<<std::setw(11)<<livetime80
		<<std::setprecision(1)<<std::setw(10)<<elivetime80<<std::endl; 

  double sumweight = 0.;
  double num = 0.;
  double sumnum = 0.;
  double avlivetime = 0.;
  double eavlivetime = 0.;
  for(unsigned int hh = 10; hh < 70; ++hh)
    {
      if(ncounter[hh] != 0)
	{
	  num = livetime[hh] * weight[hh];
	  sumnum += num;
	  sumweight += weight[hh];
	}
    }
  sumnum += livetime10*weight10 + livetime70*weight70 + livetime80*weight80;
  sumweight += weight10 + weight70 + weight80;

  avlivetime = sumnum/sumweight;
  eavlivetime = 1./(sqrt(sumweight));

  livetime_fout.unsetf(std::ios::scientific);
  livetime_fout.setf(std::ios::fixed|std::ios::floatfield);
  livetime_fout << std::endl;
  livetime_fout << std::setprecision(3) << "Average livetime = (" <<avlivetime 
		<< " pm " << std::setprecision(1) << eavlivetime << ") sec\n";
  livetime_fout << std::setprecision(3) << "                 = (" 
		<< avlivetime*sTOdays() << " pm " << std::setprecision(1)
		<< eavlivetime*sTOdays() << ") days\n";
  livetime_fout << std::endl;

  /* insert the average livetime in the header */
  fout.setf(std::ios::fixed|std::ios::floatfield);
  fout.seekp(position) << std::setprecision(3) << avlivetime << "    " 
		       << std::setprecision(1) << eavlivetime <<std::endl;

  fout.close();
  livetime_fout.close();

  return 0;
}
