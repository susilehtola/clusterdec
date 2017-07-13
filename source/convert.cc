/*
 Written by Susi Lehtola, 2016
 Copyright (c) 2016, Susi Lehtola

 This program is free software; you can redistribute it and/or
 modify it under the terms of the GNU General Public License
 as published by the Free Software Foundation; either version 2
 of the License, or (at your option) any later version.
*/

#include "parser.hh"
#include <cmath>
#include <cfloat>
#include <cstdlib>
#include <cstdio>
#include <stdexcept>
#include <algorithm>
#include <sstream>
#include <chrono>

int main(int argc, char **argv) {
  if(argc!=3) {
    printf("%s input output\n",argv[0]);
    return 1;
  }

  // Read in determinants
  size_t Norb, Nalpha, Nbeta;

  printf("Running determinant conversion \"%s\" -> \"%s\"\n",argv[1],argv[2]);
  fflush(stdout);
  
  auto tstart(std::chrono::system_clock::now());
  std::vector<determinant_string_t> detstrs(read_determinant_strings(argv[1],Norb,Nalpha,Nbeta));
  auto tstop(std::chrono::system_clock::now());

  std::chrono::duration<double> elapsed = tstop-tstart;
  printf("%i determinant strings read in %.3f seconds.\n",(int) detstrs.size(),elapsed.count());
  fflush(stdout);
  tstart=tstop;

  // Write out determinants
  write_determinant_strings(argv[2],detstrs,Norb,Nalpha,Nbeta);
  tstop=std::chrono::system_clock::now();

  elapsed = tstop-tstart;
  printf("%i determinant strings written in %.3f seconds.\n",(int) detstrs.size(),elapsed.count());
  fflush(stdout);

  return 0;
}  
