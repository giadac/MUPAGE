#ifndef DECODE_h
#define DECODE_h

#include<fstream>
#include<string>

class Decode
{
public:
  Decode(){};
  virtual ~Decode(){};

  int CommandLine(int argc, char** argv);
  int DecodeParameters();
  std::string ErrorCommandLine();
  std::string ErrorOpenFiles(int& nerror);
  std::string ErrorDecoder(int& nerror);
};

#endif /*DECODE_h*/
