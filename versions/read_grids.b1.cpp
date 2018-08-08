

#include <stdlib.h>
#include <dirent.h>
#include <stdio.h>
#include <iostream>
#include <cstring>
#include <string>
#include <vector>

#define UHBD_HEADER_LENGTH 160
#define UHBD_HEADER_TITLE_LENGTH 72
#define UHBD_BLOCK_HEADER_LENGTH 6
// calc UHBD_BLOCK_LENGTH XXXXX
#define UHBD_TAIL_LENGTH 4

/*
          June 20, 2015
    This file was created by John Ossyra for the 
    purpose of merging UHBD-format binary grids.
    These are output from the SDA BD Simulations, 
    tracking the Occupancy Probability 
    Distrubution in the Simulation Geometry. The
    workflow creates multiple executions of SDA
    so gives multiple, and equivalent PDF Grids.
    These are combined into a single averaged
    PDF. Refer to following sources:

    http://www.ccs.neu.edu/course/com3620/projects/pbe/apbs-0.2.1/tools/mesh/uhbd_asc2bin.f
    http://www.ks.uiuc.edu/Research/vmd/plugins/doxygen/uhbdplugin_8C-source.html

    compile:
       $ g++ scripts/read_grids.b1.cpp -Wall -Wno-sign-compare -o read_grids
       $ g++ scripts/read_grids.b1.cpp -Wall -Wextra -Wno-sign-compare -g -O0 -fsanitize=address -o read_grids

    run with test case:
       $ ./read_grids goodone s3s5-4adv.s3-2i2p

*/

struct UHBD_Header {
  char buffer[UHBD_HEADER_LENGTH];
};

struct UHBD_Block_Data {
  // this one determined at run time
  float* block;
};

struct UHBD_Block_Header {
  // note the kk, xm, ym are vals 2,3,4
  // rest are repeating ints pattern
  int vals[UHBD_BLOCK_HEADER_LENGTH];
};

struct UHBD_Pars {
  UHBD_Header raw;
  char  title[UHBD_HEADER_TITLE_LENGTH];
  int   xm, ym, zm;
  float ox, oy, oz;
  float scale = 1;
  float res;
};

void Write_UHBD_Header(std::string p_grid, UHBD_Pars* header) {
  FILE* f_grid;
  if ((f_grid = fopen(p_grid.c_str(), "wb")) == NULL) {
    printf("Unable to open file: %s\n", p_grid.c_str());
  }
  fwrite(header->raw.buffer, sizeof(header->raw.buffer), 1, f_grid);
  fclose(f_grid);
}

std::string Convert_Int_binString(unsigned int val) {
  std::string bin;
  unsigned int mask = 1 << (sizeof(int) * 8 - 1);
  for(int i=0; i<sizeof(int)*8; i++) {
    if((val & mask) == 0) {
      bin += '0';
    } else {
      bin += '1';
    }
    mask  >>= 1;
  }
  return bin;
}

UHBD_Pars* Open_UHBD_Read_Header(std::string p_grid) {
  FILE* f_grid;
  //deallocate sometime
  UHBD_Pars* uhbd_pars_p = new UHBD_Pars;
  printf("\nFile:  %s\n", p_grid.c_str());
  int   i_pars[8];
  float f_pars[5];

  if ((f_grid = fopen(p_grid.c_str(), "rb")) == NULL) {
    printf("Unable to open file: %s\n", p_grid.c_str());
  }
  // try to read in header, equivalent of 160 chars
  // close if unable to read in header
  if (fread(uhbd_pars_p->raw.buffer, sizeof(uhbd_pars_p->raw.buffer), 1, f_grid) != 1) {
    printf("\nIncomplete header in file:  %s\n", p_grid.c_str());
    printf("Size of header read:        %lu\n", sizeof(uhbd_pars_p->raw.buffer));
    fclose(f_grid);
  }
  // copy into respective formatted arrays
  fclose(f_grid);
  memcpy(uhbd_pars_p->title,  uhbd_pars_p->raw.buffer, sizeof(uhbd_pars_p->title));
  memcpy(i_pars, uhbd_pars_p->raw.buffer + UHBD_HEADER_TITLE_LENGTH + 8, 32);
  memcpy(f_pars, uhbd_pars_p->raw.buffer + UHBD_HEADER_TITLE_LENGTH + 40, 20);
  // then send to respective struct property
  uhbd_pars_p->xm  = i_pars[5];
  uhbd_pars_p->ym  = i_pars[6];
  uhbd_pars_p->zm  = i_pars[7];
  uhbd_pars_p->res = f_pars[1];
  uhbd_pars_p->ox  = f_pars[2];
  uhbd_pars_p->oy  = f_pars[3];
  uhbd_pars_p->oz  = f_pars[4];



//
//
//
//

  fseek(f_grid, UHBD_TAIL_LENGTH, SEEK_CUR);

  UHBD_Block_Header* block_header = new UHBD_Block_Header;

  if (fread(block_header->vals, sizeof(int), UHBD_BLOCK_HEADER_LENGTH, f_grid) != sizeof(int) * UHBD_BLOCK_HEADER_LENGTH) {
    printf("block header error, print file and block number!!\n");
  }
  printf("\nblock header values read:");
  for (int j=0; j<UHBD_BLOCK_HEADER_LENGTH; j++) {
    printf("   %d", block_header->vals[j]);
  }
  printf("\n");




//
//
//
//


  return uhbd_pars_p;
}

void Merge_UHBD_Grids (std::string d_run_home, std::string p_merged_grids, UHBD_Pars* uhbd_pars) {
  DIR* dir = opendir(d_run_home.c_str());
  std::vector<FILE*> vec_f_grids;
  struct dirent* entry = readdir(dir);
  std::vector<std::string> vec_p_grids;
  int n_entries_block;
  int n_grids;

  while (entry != NULL) {
    std::string d_name = entry->d_name;
    if (d_name.find("set") == 0){
      printf("%s\n", entry->d_name);
      vec_p_grids.push_back(d_name);
    }
    entry = readdir(dir);
  }
  n_entries_block = uhbd_pars->xm * uhbd_pars->ym;
  n_grids = vec_p_grids.size();


  // getting headers here, could do a check that all 
  // are the same, but for now just running through
  // to align for the block reads
  //
  // more importantly: building vec to store open
  // grid files, and building vec for block headers
  // though this isn't a senesible operation since
  // the format is to only accomodate the first of 
  // of these. can also ignore block headers after
  // a check that the kk-index in the file matches 
  // the z-slice for each block header
  for (int i=0; i<n_grids; i++) {
//    printf("%s, ", p_grids.at(i).c_str());
    UHBD_Header raw_header;
    std::string p_grid = d_run_home + vec_p_grids.at(i) + "/resid3d.bin.grd";
    vec_f_grids.push_back(fopen(p_grid.c_str(), "rb"));

    fread(raw_header.buffer, sizeof(raw_header.buffer), 1, vec_f_grids[i]);
    fseek(vec_f_grids[i], UHBD_TAIL_LENGTH, SEEK_CUR);
    printf("'");
    for (int j=0; j<UHBD_HEADER_TITLE_LENGTH; j++) {
      printf("%c", raw_header.buffer[j]);
    }
    printf("'\nsize of raw header buffer: %d\n", sizeof(raw_header.buffer));
  }
//// BEGIN
  //  CURRENT WORK ZONE:
  //
  // now we read in the block headers and actual
  // blocks of data. again, can throw away the 
  // headers after a check of the values.
  // 
  //             STATUS-
  //                      check values from blocks
  //                      using set123 ascii grid to navigate 'resid3d.grd'
  //                      first non-zero value:
  //                                             i = 125,230
  //                                           val = 0.47862E-03
  // 
  for (int k = 1; k <= uhbd_pars->zm; k++) {
//// current test zone
    UHBD_Block_Data* block_merged = new UHBD_Block_Data;
    block_merged->block = new float[n_entries_block]();
    for (int i=0; i<n_grids; i++) {

      UHBD_Block_Header* block_header = new UHBD_Block_Header;
      UHBD_Block_Data* block_data = new UHBD_Block_Data;
      block_data->block = new float[n_entries_block];

      if (fread(block_header->vals, sizeof(int), UHBD_BLOCK_HEADER_LENGTH, vec_f_grids[i]) != UHBD_BLOCK_HEADER_LENGTH) {
        printf("block header error, print file and block number!!\n");
      }
      printf("header values read:");
      for (int j=0; j<UHBD_BLOCK_HEADER_LENGTH; j++) {
        printf("   %d", block_header->vals[j]);
      }
      printf("\n");
      // note this is a difference in the format from
      // the sources. i find that ththis is a tail 
      // of the header, not of the block, based on 
      // ascii convresions of the grid files manually 
      // examined. SDA utilities used to convert the
      // grids show that values are 1 position right
      // meaning a value was previously read. seeking
      // before reading the block aligns the values
      // with the observable positions (checked first 
      // 2 blocks).
      fseek(vec_f_grids[i], UHBD_TAIL_LENGTH, SEEK_CUR);
      fread(block_data->block, sizeof(float), n_entries_block, vec_f_grids[i]);
      if (k == 1) {
        printf("now at file %d block %d\n", i, k);
      }
      for (int j=0; j<n_entries_block; j++) {
        if (block_data->block[j] > 1e-10) {
          block_merged->block[j] += block_data->block[j] / 128.;
          if (k == 1) {
            printf("%d:  %g, from   %g\n", j, block_merged->block[j], block_data->block[j]);
          }
        }
      }
    }
  }
  for (int i=0; i<n_grids; i++) {
    fclose(vec_f_grids[i]);
  }
  closedir(dir);
  return;
}


int main(int argc, char ** argv)
{
  // stuff used to get oriented with the set of grids
  std::string nm_group;
  std::string nm_run;
  std::string d_merge_home;
  std::string d_run_home;
  std::string d_header;
  std::string p_header;

  d_merge_home = "./";
  nm_group = argv[1];
  nm_run = argv[2];
  d_run_home = d_merge_home + "runs/" + nm_group + "/" + nm_run + "/";
  d_header = d_run_home + "set1/";
  p_header = d_header + "resid3d.bin.grd";

  // read in one header to get the parameters for each in set of grids
//  UHBD_Pars* uhbd_pars = Open_UHBD_Read_Header(p_header);
  std::string p_merged_grids;
  p_merged_grids = d_run_home + "resid3d.bin.grd.0";
  UHBD_Pars* uhbd_pars = Open_UHBD_Read_Header(p_merged_grids);

  printf("\nHeader check: z-origin= %g\n'", uhbd_pars->oz);
  for (int i=0; i<sizeof(uhbd_pars->title); i++) {
    printf("%c", uhbd_pars->title[i]);
  }
  printf("'\n");
// haven't figured out why following line doesn't work
// but all is well as shown in above loop...
//  printf("'\n'%s'\n", uhbd_pars->title);
/*
  for (int i=0; i<sizeof(uhbd_pars->raw.buffer); i++) {
    if (i < UHBD_HEADER_TITLE_LENGTH) {
      printf("buffer: %c\n", uhbd_pars->raw.buffer[i]);
    } else {
      printf("buffer: %d\n", uhbd_pars->raw.buffer[i]);
    }
  }
*/
  // write header of new uhbd file with merged grid values
//  Write_UHBD_Header(p_merged_grids, uhbd_pars);
/*
  UHBD_Pars* uhbd_pars1 = Open_UHBD_Read_Header(p_merged_grids.c_str());
  printf("\nheader check- of z-origin: %g\n", uhbd_pars1->oz);
  for (int i=0; i<sizeof(uhbd_pars1->raw.buffer); i++) {
    if (i < 72) {
      printf("buffer1: %c\n", uhbd_pars1->raw.buffer[i]);
    } else {
      printf("buffer1: %d\n", uhbd_pars1->raw.buffer[i]);
    }
  }
*/
  Merge_UHBD_Grids(d_run_home, p_merged_grids, uhbd_pars);

  delete[] uhbd_pars;
  return 0;
}
