#include <iostream> // for reading and writing files
#include <fstream>
#include <bitset>
#include <stdlib.h>
#define LINESIZE=85

using namespace std;

// Information on UHBD format, and C++ code modelled from
// following sources:
//
//     http://www.ks.uiuc.edu/Research/vmd/plugins/doxygen/uhbdplugin_8C-source.html
// Date:
//         June 13, 2015 
//
//




void* open_uhbd_read(const char* p_grid) {
  FILE* f_grid;
  char buffer[164];
  char title[72];
  int xm, ym, zm;
  float ox, oy, oz;
  float scale = 1;
  float res;
  int i_pars[8];
  float f_pars[4];

  if ((f_grid = fopen(p_grid, "rb")) == NULL) {
    printf("Unable to open file: %s\n", p_grid);
    return 1;

  // try to read in header, equivalent of 164 chars
  // close if unable to read in header
  if (fread(buffer, 1, 160, f_grid) != 160) {
    printf("Incomplete header in file: %s\n", p_grid);
    fclose(f_grid);
    return 1;
  }
  // from www.ks.uiuc.edu/Research/vmd/plugins/doxygen/uhbdplugin_8C-source.html:
  // format of header is 72 character title, followed by:
  // scale dum2 grdflag, idum2 km one km im jm km h ox oy oz
  // The first two and last four parameters are floats, the rest ints.
  //
  // first add the buffered chars to int, float data structures
  memcpy(title,  buffer, 72);
  memcpy(&scale, buffer + 72, sizeof(float));
  memcpy(i_pars, buffer + 72 + 8, 32);
  memcpy(f_pars, buffer + 72 + 40, 16);
  // then send to respective variables
  xm  = i_pars[5];
  ym  = i_pars[6];
  zm  = i_pars[7];
  res = f_pars[0];
  ox  = f_pars[1];
  oy  = f_pars[2];
  oz  = f_pars[3];

}




int main(int argc, char ** argv)
{
  int n_start =    0;
  int c_LENGTH =   500000;


  string nm_group;
  string nm_run;
  string d_merge_home;
  string d_run_home;
  d_merge_home = "./";
  nm_group = argv[1];
  nm_run = argv[2];
  d_run_home = d_merge_home + "runs/" + nm_group + "/" + nm_run + "/";
  string d_header = d_run_home + "set1/";
  string p_header = d_header + "resid3d.bin.grd";

    ifstream infile (p_header.c_str(), ios::binary | ios::in);

    char* c_buffer = new char[c_LENGTH+1];


    infile.read(c_buffer, c_LENGTH);
    c_buffer[c_LENGTH] = '\0';


    for(int i = n_start; i < c_LENGTH; i++) {
       std::bitset<8> bitform((unsigned char)(c_buffer[i]));
       std::string bitstring = bitform.to_string();
//       printf("Buffer[%d]\t is %c\n", i, (unsigned char)(c_buffer[i]));
       if ((int)(unsigned char)(c_buffer[i]) != 0) {
         printf("Buffer[%d] is\t %s\t %X\t %c\n", i, bitstring.c_str(), (unsigned char)(c_buffer[i]), (unsigned char)(c_buffer[i]));
       }
    }

    delete[] c_buffer;
    infile.close();
    return 0;
}
