#include<iostream>
#include<fstream>
#include<string>
#include<sstream>

#include "Decode.hh"
#include "Parameters.hh"

// ***********************************************************************
// Decode command line 
// ***********************************************************************
int Decode::CommandLine(int argc, char** argv)
{
  std::string sevents, sseed;
  for (int i = 1; i < argc; ++i ) /* argument 0 is the executable name */
  {
    char *arg;
    arg = argv[i];
    if (arg[0] != '-') continue;
    switch (arg[1])
      {
      case 'h' :
	std::cout << "usage: mupage.exe -[h?] <help> -i <run number> "
		  << "-n <events> -s <seed> -p <inputfile> "
		  << "-o <outputfile> <livetimeoutputfile>\n";
	exit(0);
      case 'i' :
	i++;
	numrun = std::string(argv[i]);
	break;
      case 'n' :
	i++;
	events = atoi(argv[i]);
	sevents = std::string(argv[i]);
	break;
      case 's' :
	i++;
	seed = atoi(argv[i]);
	sseed = std::string(argv[i]);
	break;
      case 'p' :
	i++;
	ParameterFileName = std::string(argv[i]);
	break;
      case 'o' :
	i++;
	OutputFileName = std::string(argv[i]);
	i++;
	LivetimeFileName = std::string(argv[i]);
	break;
      } 
  }
  if(numrun.empty() || sevents.empty() || sseed.empty() || 
     ParameterFileName.empty() || OutputFileName.empty() || 
     LivetimeFileName.empty() )
    return -79;
  else
    return 0;
}


// ***********************************************************************
// Decode parameters from the external file
// Set variables to their default values
// ***********************************************************************
int Decode::DecodeParameters()
{
  std::ifstream finpar;
  finpar.open(ParameterFileName.c_str(), std::ios::in);
  if (finpar.fail() ) 
    return -99;
  std::string token;
  char comment[256];
  while(!finpar.eof()) 
    {
      finpar >> token;
      if(token.find("#") != std::string::npos) 
	{
	  finpar.getline(comment, 255);
	  continue;
	}
      if(token == "DEPTHmax") 
	finpar >> DEPTHmax;
      else if(token == "Zmin")
	finpar >> Zmin;
      else if(token == "Zmax")
	finpar >> Zmax;
      else if(token == "CANr")
	finpar >> CANr;
      else if(token == "EnlargedCANr")
	finpar >> EnlargedCANr;
      else if(token == "THETAmin")
	{
	  finpar >> THETAmin;
	  if (THETAmin < 0.)
	    return 1;
	}
      else if(token == "THETAmax")
	{
	  finpar >> THETAmax;
	  if (THETAmax > 85.)
	    return 2;
	}
      else if(token == "Rmin")
	{
	  finpar >> Rmin;
	  if (Rmin < 0)
	    return 3;
	}
      else if(token == "Rmax")
	finpar >> Rmax;
      else if(token == "MULTmin")
	{
	  finpar >> MULTmin;
	  if (MULTmin < 1)
	    return 4;
	}
      else if(token == "MULTmax")
	finpar >> MULTmax;
      else if(token == "GEANTid")
	finpar >> GEANTid;
      else if(token == "density")
	finpar >> density;
      else if(token == "AbsLength")
	finpar >> AbsLength;
      else if(token == "NAbsLength")
	finpar >> NAbsLength;
      else if(token == "Ethreshold")
	{
	  finpar >> Ethreshold;
	  if (float(Ethreshold) < Emin)
	    return 5;
	}
    }
  finpar.close();
  return 0;
}


// ***********************************************************************
// Handling error when read command line
// ***********************************************************************  
std::string Decode::ErrorCommandLine()
{
  std::ostringstream errorcl;
  errorcl << "\nError reading command line.\n"
	  << "Check if all the options are correct typing: "
	  << "./mupage.exe -h\n";
  errorcl << std::endl;
  std::string serror = errorcl.str();
  return serror;
}

// ***********************************************************************
// Handling error when opening files
// ***********************************************************************  
std::string Decode::ErrorOpenFiles(int& nerror)
{
  std::ostringstream erroropening;
  erroropening << "\nError opening ";
  switch (nerror)
    {
    case -99:
      erroropening << "input file "<<ParameterFileName.c_str()<<std::endl;
      break;
    case -98:
      erroropening << "output file "<<OutputFileName.c_str()<<std::endl;
      break;
    case -97:
      erroropening << "output file "<<LivetimeFileName.c_str()<<std::endl;
      break;
    }
  erroropening << std::endl;
  std::string serror = erroropening.str();
  return serror;
}


// ***********************************************************************
// Handling error
// ***********************************************************************
std::string Decode::ErrorDecoder(int& nerror)
{
  std::ostringstream errordecoder;
  errordecoder << "\nERROR! \n";
  switch (nerror)
    {
    case 1:
      errordecoder << "Nonnegative input parameter for THETAmin expected.\n";
      break;
    case 2:
      errordecoder << "The chosen value of maximum zenith angle THETAmax is " 
		   << THETAmax << " degrees.\n"
		   << "The parametric formulas are valid only in the range"
		   << " of 0 to 85 degrees. \n";
      break;
    case 3:
      errordecoder << "Nonnegative input parameter for Rmin expected.\n";
      break;
    case 4:
      errordecoder << "The chosen value of minimum muon multiplicity MULTmin "
		   << "is " << MULTmin << ".\nThis value must be at least "
		   << "1.\n";
      break;
    case 5:
      errordecoder << "The chosen value of threshold energy (Ethreshold = " 
                   << Ethreshold <<" TeV) cannot be \nsmaller than the "
                   << "minimum muon energy Emin (Emin = 0.001 TeV).\n";
      break;
    case 6:
      errordecoder << "The chosen value of minimum vertical depth \n"
		   << "(i.e., the vertical depth of CAN upper disk"
		   << "with respect to the sea level) \n"
		   << "is " << DEPTHmin() << " km,"
		   << " equivalent to " << weDEPTHmin() << " km.w.e.\n"
		   << "The parametric formulas are valid only in the range"
		   << " of 1.5 to 5 km.w.e.\n"
		   << "Maybe the vertical depth of the seabed DEPTHmax \n"
		   << "is not enough deep.\n";
      break;
    case 7:
      errordecoder << "The chosen value of maximum vertical depth \n"
		   << "(seabed = DEPTHmax) is " << DEPTHmax << " km,"
		   << " equivalent to " << weDEPTHmax() << " km.w.e.\n"
		   << "The parametric formulas are valid only in the range"
		   << " of 1.5 to 5 km.w.e.\n";
      break;
    }
  errordecoder << "Change the value in file "
	       << ParameterFileName.c_str() << std::endl;
  errordecoder << std::endl;
  std::string serror = errordecoder.str();
  return serror;
}
