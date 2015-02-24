////////////Reads Unified Data Format Version 0 (UDFv0) and outputs a text file///////////// 
////////////with relevant data values///////////////////////////////////////////////////////
////////////By Tia Plautz; Modified from code written by Ford Hurley////////////////////////
////////////compile with g++ readBinary.C -o readBinary/////////////////////////////////////
////////////run with ./readBinary filename.ext//////////////////////////////////////////////


#include <stdlib.h>
#include <stdio.h>
#include <fstream>
#include <iostream>
#include <cstring>
#include <cmath>
#include <cstdarg>

using std::cout;	using std::endl;
using namespace std;


void readBinary(ifstream &infile)
{
  //Check File Type
  char magic_number[5];
  infile.read(magic_number, 4);
  magic_number[4] = '\0';
  if( strcmp(magic_number, "PCTD") ) {
    printf("Error: unknown file type (should be PCTD)!\n");
    exit(1);
  }

  int version_id;
  infile.read((char*)&version_id, sizeof(int));
  if(version_id != 0){
    printf("ERROR: Unsupported Format version (%d)!\n", version_id);
    exit(1);
  }
  
  int num_histories;
  infile.read((char*)&num_histories, sizeof(int));
  cout << "num_histories = " << num_histories << endl;
  
  printf("Reading headers from file...\n");
  
  float projection_angle, beam_energy;
  int generation_date, preprocess_date;
  int phantom_name_size, data_source_size, prepared_by_size;
  char *phantom_name, *data_source, *prepared_by;

  
  infile.read((char*)&projection_angle, sizeof(float));
  infile.read((char*)&beam_energy, sizeof(float));
  infile.read((char*)&generation_date, sizeof(int));
  infile.read((char*)&preprocess_date, sizeof(int));
  infile.read((char*)&phantom_name_size, sizeof(int));
  phantom_name = (char*)malloc(phantom_name_size);
  infile.read(phantom_name, phantom_name_size);
  infile.read((char*)&data_source_size, sizeof(int));
  data_source = (char*)malloc(data_source_size);
  infile.read(data_source, data_source_size);
  infile.read((char*)&prepared_by_size, sizeof(int));
  prepared_by = (char*)malloc(prepared_by_size);
  infile.read(prepared_by, prepared_by_size);

  cout<< "projection_angle = " << projection_angle <<endl;
  cout<< "beam_energy = " << beam_energy <<endl;
  cout<< "generation_date = " << generation_date <<endl;
  cout<< "preprocess_date = " << preprocess_date <<endl;
  cout<< "phantom_name = " << phantom_name <<endl;
  cout<< "data_source = " << data_source <<endl;
  cout<< "prepared_by = " << prepared_by <<endl;
  
  printf("Loading %d histories from file\n", num_histories);
  
  int data_size = num_histories * sizeof(float);

  /////////////////////

  int max_hist_per_projection = num_histories;
  unsigned int mem_size_hist_floats = sizeof(float) * max_hist_per_projection;
  unsigned int mem_size_hist_ints = sizeof(int) * max_hist_per_projection;
  
  // Allocate memory
  float* t_in1_h         = (float*) malloc(mem_size_hist_floats);
  float* t_in2_h         = (float*) malloc(mem_size_hist_floats);
  float* t_out1_h        = (float*) malloc(mem_size_hist_floats);
  float* t_out2_h        = (float*) malloc(mem_size_hist_floats);
  float* u_in1_h         = (float*) malloc(mem_size_hist_floats);
  float* u_in2_h         = (float*) malloc(mem_size_hist_floats);
  float* u_out1_h        = (float*) malloc(mem_size_hist_floats);
  float* u_out2_h        = (float*) malloc(mem_size_hist_floats);
  float* v_in1_h         = (float*) malloc(mem_size_hist_floats);
  float* v_in2_h         = (float*) malloc(mem_size_hist_floats);
  float* v_out1_h        = (float*) malloc(mem_size_hist_floats);
  float* v_out2_h        = (float*) malloc(mem_size_hist_floats);
  float* WEP_h           = (float*) malloc(mem_size_hist_floats);
  int* lucy_rot_h = (int*) malloc(mem_size_hist_ints);

  
  // Initialize memory
  memset(t_in1_h, 0, mem_size_hist_floats);
  memset(t_in2_h, 0, mem_size_hist_floats);
  memset(t_out1_h, 0, mem_size_hist_floats);
  memset(t_out2_h, 0, mem_size_hist_floats);
  memset(v_in1_h, 0, mem_size_hist_floats);
  memset(v_in2_h, 0, mem_size_hist_floats);
  memset(v_out1_h, 0, mem_size_hist_floats);
  memset(v_out2_h, 0, mem_size_hist_floats);
  memset(u_in1_h, 0, mem_size_hist_floats);
  memset(u_in2_h, 0, mem_size_hist_floats);
  memset(u_out1_h, 0, mem_size_hist_floats);
  memset(u_out2_h, 0, mem_size_hist_floats);
  memset(WEP_h, 0, mem_size_hist_floats);
  memset(lucy_rot_h, 0, mem_size_hist_ints);
  /////////////////////
  
  infile.read((char*)t_in1_h, data_size);
  infile.read((char*)t_in2_h, data_size);
  infile.read((char*)t_out1_h, data_size);
  infile.read((char*)t_out2_h, data_size);
  infile.read((char*)v_in1_h, data_size);
  infile.read((char*)v_in2_h, data_size);
  infile.read((char*)v_out1_h, data_size);
  infile.read((char*)v_out2_h, data_size);
  infile.read((char*)u_in1_h, data_size);
  infile.read((char*)u_in2_h, data_size);
  infile.read((char*)u_out1_h, data_size);
  infile.read((char*)u_out2_h, data_size);
  infile.read((char*)WEP_h, data_size);
  
  
  // for( int i = 0; i < num_histories; i++ )
  for( int i = 0; i<num_histories && i<10; i++ )
  {
    // convert to cm:
    t_in1_h[i]  /= 10.0;
    t_in2_h[i]  /= 10.0;
    t_out1_h[i] /= 10.0;
    t_out2_h[i] /= 10.0;
    v_in1_h[i]  /= 10.0;
    v_in2_h[i]  /= 10.0;
    v_out1_h[i] /= 10.0;
    v_out2_h[i] /= 10.0;
    u_in1_h[i]  /= 10.0;
    u_in2_h[i]  /= 10.0;
    u_out1_h[i] /= 10.0;
    u_out2_h[i] /= 10.0;
    WEP_h[i] /= 10.0;
    lucy_rot_h[i] = projection_angle;


    
    printf("%8.4f \t %8.4f \t %8.4f \t %8.4f \t %8.4f \t %8.4f \t %8.4f \t %8.4f \t %8.4f \n",   //
	   v_in1_h[i], v_in2_h[i],  v_out1_h[i], v_out2_h[i],  t_in1_h[i],  t_in2_h[i], //
	   t_out1_h[i], t_out2_h[i], WEP_h[i]);		

  }
  
}


// int main(){
  
//   const char directory[] = "/Users/Tia/pct-sim";
//   char ifname[512];
//   sprintf(ifname, "%s/pCTraw_Run_51.out-0.root.bin", directory);
//   ifstream projFile(ifname, ios::binary);
//   readBinary(projFile);

//   return 0;

// }


int main(int argc, char *argv[])
{
  //cout<< "argc = " << argc <<endl;  for (int i=0; i<argc; ++i) cout<< i <<" "<< argv[i] <<endl;

  if (argc == 1) {
    cout<< "Usage:\n" << argv[0] << " par1 par2 par3" <<endl;
    return 0;
  }

  const char* fileName = argv[1];

  ifstream projFile(fileName, ios::binary);
  readBinary(projFile);

  return 0;
}
