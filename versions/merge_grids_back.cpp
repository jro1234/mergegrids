
#include <unistd.h>


#include <stdlib.h>
#include <dirent.h>
#include <stdio.h>
#include <omp.h>
#include <iostream>
#include <cstring>
#include <string>
#include <vector>

#define UHBD_HEADER_LENGTH 160
#define UHBD_HEADER_TITLE_LENGTH 72
// the final definition of block header lenghts
// are system dependent on the system that wrote the original file.
// make sure this is integrated to run on the 
// original system, immediately consolidating 
// these output grid files. 
// luckily most have sizeof(int) === 4.
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
       $ g++ scripts/merge_grids.cpp -fopenmp -std=c++11 -Wall -Wno-sign-compare -o merge_grids
       $ g++ scripts/merge_grids.cpp -fopenmp -std=c++11 -Wall -Wextra -Wno-sign-compare -g -O0 -fsanitize=address -o merge_grids

    run with test case:
       $ ./merge_grids goodone s3s5-4adv.s3-2i2p

*/


// some on internet say gnu c++ compiler
// doesn't work with omp_get_num_threads
// here is replacement, detecting threads
// by sinking each into reduction of 'n':
int omp_get_n_threads() {
  int n = 0;

  #pragma omp parallel reduction(+:n)
  n += 1;
  return n;
}


struct UHBD_Header {
  char raw[UHBD_HEADER_LENGTH + UHBD_TAIL_LENGTH];
  int  xm, ym, zm;
//  char  title[UHBD_HEADER_TITLE_LENGTH];
  float ox;
  float oy;
  float oz;
//  float scale = 1;
//  float res;
};

struct UHBD_Block {
  float  factor = {0};
// extra calculations
//  float  sum = {0};
  int    n_entries;
  int    kk;
//  int    xm;
//  int    ym;
  float* data;
  int    header[UHBD_BLOCK_HEADER_LENGTH];
  char   raw_header[sizeof(int) * UHBD_BLOCK_HEADER_LENGTH + UHBD_TAIL_LENGTH];
  int    id = -1;
};

struct UHBD_Grid {
  UHBD_Header header;
  std::vector<UHBD_Block*> blocks;
  int   n_entries_block;
//  float sum = {0};
  int   n_blocks;
  float factor = {0};
  int   id = -1;
};


void Print_UHBD_Header (UHBD_Header header) {
  float firstf[2];
  float secondf[10];
  int firsti[8];
  int secondi[5];

  int offset;

  printf("\nPrinting out the header\n");
  printf("'");
  for (int i=0; i<UHBD_HEADER_TITLE_LENGTH; i++) {
    printf("%c", header.raw[i]);
  }
  printf("'");

  offset = UHBD_HEADER_TITLE_LENGTH;
  memcpy(firstf, header.raw + offset, 2*sizeof(float));
  printf("\n%f", firstf[0]);
  printf("\n%f", firstf[1]);

  offset += 2*sizeof(float);
  memcpy(firsti, header.raw + offset, 8*sizeof(int));
  for (int i=0; i<8; i++) {
    printf("\nAt byte #%u", offset + i * sizeof(int));
    printf("\n%d", firsti[i]);
  }

  offset += 8*sizeof(int);
  memcpy(secondf, header.raw + offset, 10*sizeof(float));
  for (int i=0; i<10; i++) {
    printf("\n%f", secondf[i]);
  }

  offset += 10*sizeof(float);
  memcpy(secondi, header.raw + offset, 2*sizeof(int));
  for (int i=0; i<2; i++) {
    printf("\n%d", secondi[i]);
  }

  offset += 2*sizeof(int);
  printf("\nTotal bytes read so far: %d\n", offset);
}


UHBD_Grid* Read_UHBD_Header (FILE* f_grid) {
  UHBD_Grid* data_grid = new UHBD_Grid;
  int   i_pars[8];
  float f_pars[3];

  // try to read in header, equivalent
  // of 160 chars. close if unable to 
  // read in header. 
  // there is a tail to header as well.
  if (fread(data_grid->header.raw, UHBD_HEADER_LENGTH + UHBD_TAIL_LENGTH, 1, f_grid) != 1) {
    printf("\nIncomplete header in file:  \n");
    printf("Size of header read:        %d\n", UHBD_HEADER_LENGTH + UHBD_TAIL_LENGTH);
    fclose(f_grid);
  } 
  else {
    memcpy(i_pars, data_grid->header.raw + UHBD_HEADER_TITLE_LENGTH + 8, 32);
    memcpy(f_pars, data_grid->header.raw + UHBD_HEADER_TITLE_LENGTH + 4*sizeof(float) + 8*sizeof(int), 3*sizeof(float));

    // then send to respective struct property
    data_grid->header.xm = i_pars[5];
    data_grid->header.ym = i_pars[6];
    data_grid->header.zm = i_pars[7];
    data_grid->header.ox = f_pars[0];
    data_grid->header.oy = f_pars[1];
    data_grid->header.oz = f_pars[2];
    data_grid->n_entries_block = data_grid->header.xm * data_grid->header.ym;
    data_grid->n_blocks = data_grid->header.zm;
  }

  return data_grid;
}


UHBD_Grid* Initialize_Target_Grid (UHBD_Grid* data_grid, int resize_grid) {
  UHBD_Grid* target_grid = new UHBD_Grid;
  int i_pars[8];

  memcpy(target_grid->header.raw, data_grid->header.raw, UHBD_HEADER_LENGTH + UHBD_TAIL_LENGTH);

  if (resize_grid > 0) {
    int first_chunk_size;
    int offset;
    int nx;
    first_chunk_size = UHBD_HEADER_TITLE_LENGTH + 2*sizeof(float) + 3*sizeof(int);

    offset = first_chunk_size;
    memcpy(&target_grid->header.raw[offset], reinterpret_cast<char*>(&resize_grid), sizeof(int));

    Print_UHBD_Header(target_grid->header);

    offset += sizeof(int);
    memcpy(&target_grid->header.raw[first_chunk_size+2*sizeof(int)], reinterpret_cast<char*>(&resize_grid), sizeof(int));
    memcpy(&target_grid->header.raw[first_chunk_size+3*sizeof(int)], reinterpret_cast<char*>(&resize_grid), sizeof(int));
    memcpy(&target_grid->header.raw[first_chunk_size+4*sizeof(int)], reinterpret_cast<char*>(&resize_grid), sizeof(int));

    nx = (data_grid->header.xm - resize_grid) / 2;
    printf("\n nX is:  %d", nx);
    target_grid->header.ox = data_grid->header.ox + (float)nx;
    target_grid->header.oy = data_grid->header.oy + (float)nx;
    target_grid->header.oz = data_grid->header.oz + (float)nx;
    target_grid->header.xm = resize_grid;
    target_grid->header.ym = resize_grid;
    target_grid->header.zm = resize_grid;
    target_grid->n_blocks = target_grid->header.zm;
    target_grid->n_entries_block = target_grid->header.xm * target_grid->header.ym;
    printf("\n\noXd:  %f", data_grid->header.ox);
    printf("\n\noXt:  %f", target_grid->header.ox);
    printf("\n\noYd:  %f", data_grid->header.oy);
    printf("\n\noYt:  %f", target_grid->header.oy);
    printf("\n\noZd:  %f", data_grid->header.oz);
    printf("\n\noZt:  %f", target_grid->header.oz);

    memcpy(&target_grid->header.raw[first_chunk_size+4*sizeof(int)+3*sizeof(float)], reinterpret_cast<char*>(&target_grid->header.ox), sizeof(float));
    memcpy(&target_grid->header.raw[first_chunk_size+4*sizeof(int)+4*sizeof(float)], reinterpret_cast<char*>(&target_grid->header.oy), sizeof(float));
    memcpy(&target_grid->header.raw[first_chunk_size+4*sizeof(int)+5*sizeof(float)], reinterpret_cast<char*>(&target_grid->header.oz), sizeof(float));

  }
  else {
    memcpy(i_pars, target_grid->header.raw + UHBD_HEADER_TITLE_LENGTH + 8, 32);

    target_grid->header.xm  = i_pars[5];
    target_grid->header.ym  = i_pars[6];
    target_grid->header.zm  = i_pars[7];
    target_grid->n_blocks = target_grid->header.zm;
    target_grid->n_entries_block = target_grid->header.xm * target_grid->header.ym;
  }

  ///////
  ///////
  ///////   Fixing Chimera's problem with the format. Need a value of 1 in
  ///////   this particular position for the reader to believe its UHBD grid.
  ///////
  int a_one;
  a_one = 1;
  memcpy(&target_grid->header.raw[96], reinterpret_cast<char*>(&a_one), sizeof(int));
  ///////
  ///////
  
  printf("\nPrinting target grid header, should have n 1 n n n 0.000 1.000 ox oy oz");
  Print_UHBD_Header(target_grid->header);

  for (int kk=0; kk<target_grid->n_blocks; kk++) {
    UHBD_Block* block = new UHBD_Block;
    float* data = new float[target_grid->n_entries_block]();

    // note difference between loop / file kk
    block->kk = kk + 1;
    block->data = data;
    block->n_entries = target_grid->n_entries_block;

    target_grid->blocks.push_back(block);
  }

  return target_grid;
}


void Write_UHBD_Block_Header (FILE* f_grid, UHBD_Block* block) {
  printf("\nblock header values to write:");
  for (int j=0; j<UHBD_BLOCK_HEADER_LENGTH; j++) {
    printf("   %d", block->header[j]);
  }

  printf("\n");
  fwrite(block->raw_header, sizeof(int) * UHBD_BLOCK_HEADER_LENGTH + UHBD_TAIL_LENGTH, 1, f_grid);
}


void Write_UHBD_Header (FILE* f_grid, UHBD_Header header) {
  fwrite(header.raw, UHBD_HEADER_LENGTH + UHBD_TAIL_LENGTH, 1, f_grid);
}


void Read_UHBD_Block (UHBD_Block* block, FILE* f_grid) {
    float* data = new float[block->n_entries];

//    printf("\nnew block, size: %d\n", block->n_entries);
    if (fread(block->raw_header, sizeof(int) * UHBD_BLOCK_HEADER_LENGTH + UHBD_TAIL_LENGTH, 1, f_grid) != 1) {
      printf("\nBlock header read error\n");
    }
    memcpy(block->header, block->raw_header, sizeof(int) * UHBD_BLOCK_HEADER_LENGTH);
    block->kk = block->header[2];

    // ???????
    // This check needs some form of change...
    // ???????
    if (block->header[2]%100==1)
    {
      printf("\nblock header values read:");
      for (int j=0; j<UHBD_BLOCK_HEADER_LENGTH; j++) {
        printf("   %d", block->header[j]);
      }
      printf("\n");
    }
    // ???????
    // ???????

    if (fread(data, sizeof(float), block->n_entries, f_grid) != block->n_entries) {
      printf("\nBlock read error\n");
    }
    else {
      block->data = data;
    }
}


void Read_UHBD_Grid (UHBD_Grid* data_grid, FILE* f_grid) {
  for (int kk=0; kk<data_grid->n_blocks; kk++) {
    UHBD_Block* block = new UHBD_Block;

    block->id = data_grid->id;
    block->factor = data_grid->factor;
    block->n_entries = data_grid->n_entries_block;
    ///////////////////////////////////////
    // development convenience structure //
    //if (kk == 0)                       //
    //{                                  //
    ///////////////////////////////////////
    Read_UHBD_Block(block, f_grid);
    ///////////////////////////////////////
    //}                                  //
    // development convenience structure //
    ///////////////////////////////////////
    data_grid->blocks.push_back(block);
  }
}


void Reduce_UHBD_Block (UHBD_Block* target_block, UHBD_Block* data_block, float reduce_factor) {
//  printf("\nreduce factor: %f", reduce_factor);
//  printf("\nn entries    : %d\n", target_block->n_entries);

  for (int i=0; i<target_block->n_entries; i++) {
    if (data_block->data[i] > 1e-10) {
      target_block->data[i] += reduce_factor * data_block->data[i];
      if (target_block->kk == 1)
      {
        printf("\n non-zero value: block %d, grid %d:  %d  %.9f\n", data_block->kk, data_block->id, i, data_block->data[i]);
        printf("   new   value : block %d, grid %d:  %d  %.9f\n", target_block->kk, target_block->id, i, target_block->data[i]);

        float old_value = target_block->data[i] - reduce_factor * data_block->data[i];
        if (old_value > 1e-10){
          printf("   old   value : block %d, grid %d:  %d  %.9f\n", target_block->kk, target_block->id, i, old_value);
        }
      }
    }
  }

  delete[] data_block->data;
  delete data_block;
}


//void Reduce_Write_UHBD_Block (FILE* f_grid, UHBD_Block* first_block, UHBD_Block* second_block, float reduce_factor) {
////  printf("\nreduce factor: %f", reduce_factor);
////  printf("\nn entries    : %d\n", target_block->n_entries);
//  Write_UHBD_Block_Header(f_grid, first_block);
//
//  delete[] second_block->data;
//  delete[] second_block;
//  delete[] first_block->data;
//  delete[] first_block;
//}


void Write_UHBD_Block (FILE* f_grid, UHBD_Block* block) {
//  printf("\nreduce factor: %f", reduce_factor);
//  printf("\nn entries    : %d\n", target_block->n_entries);
  Write_UHBD_Block_Header(f_grid, block);
  fwrite(block->data, sizeof(float) * block->n_entries, 1, f_grid);

  delete[] block->data;
  delete block;
}


void Reduce_UHBD_Grids (UHBD_Grid* target_grid, UHBD_Grid* data_grid) {
  for (int kk=0; kk<target_grid->n_blocks; kk++) {

    ///////////////////////////////////////
    ///////////////////////////////////////
    // development convenience structure //
    //if (kk == 0)                       //
    //{                                  //
    ///////////////////////////////////////
    ///////////////////////////////////////

      //????????????????????
      // See if this if condition can go outside of loop and avoid
      // checking over and over again, just set one time
      if (target_grid->blocks[kk]->factor < 1e-10) {
        memcpy(target_grid->blocks[kk]->raw_header, data_grid->blocks[kk]->raw_header, sizeof(int) * UHBD_BLOCK_HEADER_LENGTH + UHBD_TAIL_LENGTH);
        memcpy(target_grid->blocks[kk]->header, target_grid->blocks[kk]->raw_header, sizeof(int) * UHBD_BLOCK_HEADER_LENGTH);
      }
      //????????????????????

      Reduce_UHBD_Block(target_grid->blocks[kk], data_grid->blocks[kk], data_grid->blocks[kk]->factor);
      target_grid->blocks[kk]->factor += data_grid->factor;

    ///////////////////////////////////////
    ///////////////////////////////////////
    //}                                  //
    // development convenience structure //
    ///////////////////////////////////////
    ///////////////////////////////////////

  }

  target_grid->factor += data_grid->factor;
}

void Merge_UHBD_Block (std::vector<UHBD_Grid*> vec_grids, int n_block, int n_grids) {
  vec_grids[0]->id = 0;
  for (int i=1; i<n_grids; i++) {
    vec_grids[i]->id = i;
    // reduce to 0th to store new values.
    // note all values are normalized to
    // their final contribution to the 
    // averaged grid with the first 
    // reduction factor, so it's one 
    // from here on out.
    Reduce_UHBD_Block(vec_grids[0]->blocks[n_block], vec_grids[i]->blocks[n_block], 1);
  }

  vec_grids[n_grids-1]->id = n_grids - 1;
}


void Merge_UHBD_Grid (std::vector<UHBD_Grid*> vec_grids, int* parts_target_grid, int n_thread, int n_threads) {
  int blocks_interval[2];

  if (n_thread == 0) {
    blocks_interval[0] = 0;
    blocks_interval[1] = parts_target_grid[n_thread] - 1;
  }
  else {
    blocks_interval[0] = parts_target_grid[n_thread-1];
    blocks_interval[1] = parts_target_grid[n_thread] - 1;
  }

  printf("\n\n%d Now merging final grid", n_thread);
  printf("\n%d partition of blocks: [ %d,  %d ]\n", n_thread, blocks_interval[0], blocks_interval[1]);

  // since we calculated a partition of this
  // grid structure for each thread, there 
  // is no overlap in the grid-block 
  // operations in this loop.
  for (int i=blocks_interval[0]; i<=blocks_interval[1]; i++) {
    // we use the first grid in the vector
    // of 1st round merger to store the 
    // 2nd round merger. the final pass
    // will directly write to file.
    Merge_UHBD_Block(vec_grids, i, n_threads);
  }
}

void Write_UHBD_Grid (FILE* f_target_grid, UHBD_Grid* grid) {
  Write_UHBD_Header(f_target_grid, grid->header);

  for (int i=0; i<grid->n_blocks; i++) {
    Write_UHBD_Block(f_target_grid, grid->blocks[i]);
  }
}


void Merge_UHBD_Grids (int** lists_p_grids, std::vector<std::string>& vec_p_grids, std::string d_run_home, int n_threads, int n_grids_thread, int n_grids, int resize_grid) {
  std::vector<UHBD_Grid*> vec_merged_grids(n_threads);
  int parts_merged_grid[n_threads];

  omp_set_num_threads(n_threads);
  #pragma omp parallel
  {
    int n_thread = omp_get_thread_num();
    printf("launching number:  %d\n", n_thread);

    for (int i=0; i<n_grids_thread; i++) {
      int idx_p_grid = lists_p_grids[n_thread][i] - 1;

      // now done with grids for this thread
      if (idx_p_grid == -1) {
        printf("done with grids in thread number:  %d\n", n_thread);
        break;
      }

      else {
        FILE* f_grid;
        std::string p_grid = vec_p_grids[idx_p_grid];

        if ((f_grid = fopen(p_grid.c_str(), "rb")) == NULL) {
          printf("Unable to open file: %s\n", p_grid.c_str());
        }

        UHBD_Grid* data_grid = Read_UHBD_Header(f_grid);

        if (n_thread == 0) {
          Print_UHBD_Header(data_grid->header);
        }

        // initialize the data_grid which will be used
        // only locally as it is reduced into the 
        // target grid, merged_grid
        data_grid->factor = 1 / (float)n_grids;
        data_grid->id = idx_p_grid;

        printf("\n Checking data file from thread #%d", n_thread);
        printf("\n Right now, i=%d", i);
        printf("\n%d Header check: x-dim = %d", n_thread, data_grid->header.xm);
        printf("\n%d Header check: y-dim = %d", n_thread, data_grid->header.ym);
        printf("\n%d Header check: z-dim = %d", n_thread, data_grid->header.zm);
        printf("\n%d # of entries / block= %d", n_thread, data_grid->n_entries_block);
        printf("\n%d Factor check: factor= %f\n", n_thread, data_grid->factor);

        if (i == 0) {
          UHBD_Grid* merged_grid;

          // this is where the target grid is initialized
          // it is locally accessed as merged_grid.
          // each thread uses n_thread as address within
          // vec_merged_grids so merger can level up.
          merged_grid = Initialize_Target_Grid(data_grid, resize_grid);
          vec_merged_grids[n_thread] = merged_grid;

          printf("\n Checking target grid from thread #%d", n_thread);
          printf("\n The grid is in memory and should be initialized to all zeroes");
          printf("\n%d Header check: x-dim = %d", n_thread, merged_grid->header.xm);
          printf("\n%d Header check: y-dim = %d", n_thread, merged_grid->header.ym);
          printf("\n%d Header check: z-dim = %d", n_thread, merged_grid->header.zm);
          printf("\n%d # of entries / block= %d", n_thread, merged_grid->n_entries_block);
          printf("\n%d Factor check: factor= %f", n_thread, merged_grid->factor);
          printf("\n%d           :[middle] = %f\n", n_thread, merged_grid->blocks[50]->data[31615]); 
          printf("\n%d             :[last] = %f\n", n_thread, merged_grid->blocks[50]->data[249999]); 
        }

        Read_UHBD_Grid(data_grid, f_grid);

        printf("\n Checking data grid after reading from thread #%d", n_thread);
        printf("\n The data is in memory and should have some non-zero values\n");
        printf("%d     From   file : %s\n", n_thread, p_grid.c_str());
        printf("%d  data grid check: zm     =  %d\n", n_thread, data_grid->header.zm); 
        printf("%d                 : factor =  %f\n", n_thread, data_grid->factor);
        printf("%d                 : [non-0]=  %f\n", n_thread, data_grid->blocks[0]->data[127743]); 

        Reduce_UHBD_Grids(vec_merged_grids[n_thread], data_grid);

        printf("\n Checking target grid after reducing data from thread #%d", n_thread);
        printf("\n The data is in memory and should have some non-zero values\n");
        printf("%d     From   file : %s\n", n_thread, p_grid.c_str());
        printf("%dtarget grid check: zm     =  %d\n", n_thread, vec_merged_grids[n_thread]->header.zm); 
        printf("%d                 : factor =  %f\n", n_thread, vec_merged_grids[n_thread]->factor);
        //printf("%d                 : [non-0]=  %f\n", n_thread, vec_merged_grids[n_thread]->blocks[0]->data[127743]); 

        delete data_grid;
      }
    }

    printf("\n%d This thread is done reading grid data files", n_thread);

    // calculate interval for each thread
    // this algorithm creates a balanced 
    // partition of the grid blocks when 
    // applied over all threads (haven't
    // considered n_threads==1)
    // MUST USE vec_merged_grids[n_thread] to get n_blocks because the other
    //          grids in vec_merged_grids are at unknown stage of construction
    //          -> we know the merged_grid of n_thread is all the way done
    #pragma omp single
    {
      printf("\n%d 2nd stage of merge: \n\tdetermining grid partitioning\n", n_thread);

      for (int i=0; i<n_threads; i++) {
        parts_merged_grid[i] = (i + 1) * (vec_merged_grids[n_thread]->n_blocks / n_threads);

        if (i < vec_merged_grids[n_thread]->n_blocks % n_threads) {
          parts_merged_grid[i] += i + 1;
        }
        else {
          parts_merged_grid[i] += vec_merged_grids[n_thread]->n_blocks % n_threads;
        }

        printf("\n%d Merged grid parts: %d, %d", n_thread, i, parts_merged_grid[i]);
      }
    } //endof omp single

    printf("\n%d Entering final stage, all threads should be done reading grid data files", n_thread);

    Merge_UHBD_Grid(vec_merged_grids, parts_merged_grid, n_thread, n_threads);

  } //endof omp parallel

  // now create intermediate segments of final file
  // note this is not currently geared for 1 thread
  std::string p_merged_grid;
  FILE* f_merged_grid;

  p_merged_grid = d_run_home + "resid3d.merged.bin.grd";
  if ((f_merged_grid = fopen(p_merged_grid.c_str(), "wb")) == NULL) {
    printf("Unable to open file: %s\n", p_merged_grid.c_str());
  }

  Write_UHBD_Grid(f_merged_grid, vec_merged_grids[0]);
}


int Read_d_run (std::vector<std::string>& vec_p_grids, std::string d_run_home) {
  DIR* dir = opendir(d_run_home.c_str());
  struct dirent* entry = readdir(dir);

  while (entry != NULL) {
    std::string d_name = entry->d_name;
    if (d_name.find("set") == 0){
      std::string p_grid;
//      printf("%s\n", entry->d_name);
      p_grid = d_run_home + d_name + "/resid3d.bin.grd";
      vec_p_grids.push_back(p_grid);
    }
    entry = readdir(dir);
  }
  return vec_p_grids.size();
}


int main(int argc, char ** argv)
{
  // stuff used to get oriented with the set of grids
  std::string nm_group;
  std::string nm_run;
  std::string d_merge_home;
  std::string d_run_home;
  //std::string d_header;
  //std::string p_header;

  int resize_grid;

  d_merge_home = "./";
  printf("\nMerging Grids Program\n\n\n");
  printf("Building the working directory of merger: '%s'", d_merge_home.c_str());

  nm_group = argv[1];
  nm_run = argv[2];
  if (argc==3) {
    resize_grid = 0;
  } else if (argc==4) {
    resize_grid = atoi(argv[3]);
  }

  d_run_home = d_merge_home + "runs/" + nm_group + "/" + nm_run + "/";
  printf("\nBuilding the working directory of merger: '%s'\n", d_run_home.c_str());

  std::vector<std::string> vec_p_grids;
  int n_grids;
  int n_threads;
  int n_grids_thread;
  int n_grids_extra;
  int** lists_p_grids;

  n_grids = Read_d_run(vec_p_grids, d_run_home);
  n_threads = omp_get_n_threads();
  // incase of tortoises
  // also good for checks
  //n_threads = 1;

  n_grids_thread = n_grids / n_threads;
  n_grids_extra  = n_grids % n_threads;
  if (n_grids_extra > 0) {
    n_grids_thread += 1;
  }

  lists_p_grids = new int *[n_threads];
  for (int i=0; i<n_grids_thread; i++) {
    for (int j=0; j<n_threads; j++) {
      if (i == 0) {
        // format is to start at 1, then if
        // a grid has the zero index (and is
        // the final entry as a check) we 
        // know it was unfilled in the extra
        // line: below is zero-init method
        lists_p_grids[j] = new int[n_grids_thread]();
      }
      int val = i * n_threads + j + 1;

      if (val > n_grids) {
        // leave this and remaining as zero
        break;
      }
      else {
        lists_p_grids[j][i] = val;
      }
    }
  }

  printf("\n\nMerging Grids in folder:  %s\n\n\nmerge config: \n\t\tn_grids: %d    n_threads: %d    n_grids_thread: %d\n\n", d_run_home.c_str(), n_grids, n_threads, n_grids_thread);

  for (int i=0; i<vec_p_grids.size(); i++) {
    printf("  grid# %d\t%s\n", i, vec_p_grids[i].c_str());
  }

  printf("\ngrids-threads breakdown:\n");
  for (int i=0; i<n_threads; i++) {
    printf("\t\tthread# %d)   ", i);
    for (int j=0; j<n_grids_thread; j++) {
      if (lists_p_grids[i][j] > 0) {
        printf("%d  ", lists_p_grids[i][j]);
      }
    }
    printf("\n");
  }
  printf("\n");

  Merge_UHBD_Grids(lists_p_grids, vec_p_grids, d_run_home, n_threads, n_grids_thread, n_grids, resize_grid);

  printf("\n\n");

  return 0;
}
