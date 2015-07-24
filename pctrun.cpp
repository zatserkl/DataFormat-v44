/*
Linux:
g++ -Wall -o pctrun pctrun.cpp
Mac OS X:
clang++ -Wall -o pctrun pctrun.cpp
*/
#include <iostream>
#include <fstream>
#include <cstdarg>
#include <iomanip>

using std::cout;      using std::endl;

int main(int argc, char *argv[])
{
  //cout<< "argc = " << argc <<endl;  for (int i=0; i<argc; ++i) cout<< i <<" "<< argv[i] <<endl;

  if (argc == 1) {
    cout<< "Usage:\n" << argv[0] << " <name of the pCT binary file>" <<endl;
    return 0;
  }

  const char* ifname = argv[1];
  std::ifstream ifile(ifname, std::ios::binary);
  if (!ifile) {
    cout<< "Could not open file " << ifname <<endl;
    return 0;
  }
  else cout<< "\n--> check file " << ifname <<endl;

  char byte = 0;
  // for (int ibyte=0; ibyte<4; ++ibyte) {
  //   ifile.get(byte);
  //   cout<< std::setw(2) << std::setfill('0') << std::hex << (unsigned) (unsigned char) byte << std::dec << " ";
  // }
  // cout<<endl;

  const unsigned char pattern_run_header_identifier[3] = {0xD2, 0x55, 0x4E};     // 1R U N

  // run header identifier (including the first zero byte)

  unsigned char buf_run_header_identifier[4];
  for (int ibyte=0; ibyte<4; ++ibyte) {
    ifile.get(byte);
    cout<< std::setw(2) << std::setfill('0') << std::hex << (unsigned) (unsigned char) byte << std::dec << " ";
    buf_run_header_identifier[ibyte] = byte;
  }
  cout<<endl;

  bool good_run_header = true
    && buf_run_header_identifier[0] == 0
    && buf_run_header_identifier[1] == pattern_run_header_identifier[0]
    && buf_run_header_identifier[2] == pattern_run_header_identifier[1]
    && buf_run_header_identifier[3] == pattern_run_header_identifier[2]
  ;
  if (!good_run_header) {
    cout<< "This is not a pCT binary data file" <<endl;
    return 0;
  }

  union {
    char buf[4];
    int run;
  } run_union;
  run_union.run = 0;

  for (int ibyte=0; ibyte<3; ++ibyte) {
    ifile.get(byte);
    cout<< std::setw(2) << std::setfill('0') << std::hex << (unsigned) (unsigned char) byte << std::dec << " ";
    run_union.buf[2 - ibyte] = byte;
  }
  cout<<endl;
  cout<< "run = " << run_union.run <<endl;

  return 0;
}
