/* THIS FILE IS AUTOGENERATED, DON'T TOUCH IT! */

#include "gather.hh"
#include "parser.hh"
#include "tensor.hh"
#include "types.hh"
#include <vector>
#include <climits>
#include <cmath>
#include <chrono>
#include <sstream>
#include <cstdio>

/*
 Written by Susi Lehtola, 2016
 Copyright (c) 2016, Susi Lehtola

 This program is free software; you can redistribute it and/or
 modify it under the terms of the GNU General Public License
 as published by the Free Software Foundation; either version 2
 of the License, or (at your option) any later version.
*/
/* This code is autogenerated by generate.cc. Don't edit it by hand. */
int main(int argc, char **argv) {
  printf("CLUSTERDEC, written by Susi Lehtola (2016)\n\n");
  if(argc>4) {
    printf("Usage: %s (file) (maxr) (ndets)\n",argv[0]);
    return 1;
  }
  auto tstart(std::chrono::system_clock::now());
  std::string fname = (argc>1) ? argv[1] : "wfs.dat";
  int maxr = (argc>2) ? atoi(argv[2]) : 6;
  int ndets = 0;
  double thr = 1e-8;
  if(argc>3) {
    if(atof(argv[3])>1) {
      ndets=atoi(argv[3]);
    } else {
      thr=atof(argv[3]);
    }
  }
  std::vector<determinant_t> list(parse(fname,ndets,thr));
  std::vector<double> Cnorm(17,0.0), Tnorm(9,0.0);

  // Form excitation tensors
  Tensor1 T1(list);
  Tensor2 T2(list);
  Tensor3 T3(list);
  Tensor4 T4(list);
  Tensor5 T5(list);
  Tensor6 T6(list);
  Tensor7 T7(list);
  Tensor8 T8(list);
  Tensor9 T9(list);
  Tensor10 T10(list);
  Tensor11 T11(list);
  Tensor12 T12(list);
  Tensor13 T13(list);
  Tensor14 T14(list);
  Tensor15 T15(list);
  Tensor16 T16(list);

  // Check maximum rank
  int nzt=0;
  if(T8.size()) {
    nzt=8;
  } else if(T7.size()) {
    nzt=7;
  } else if(T6.size()) {
    nzt=6;
  } else if(T5.size()) {
    nzt=5;
  } else if(T4.size()) {
    nzt=4;
  } else if(T3.size()) {
    nzt=3;
  } else if(T2.size()) {
    nzt=2;
  } else if(T1.size()) {
    nzt=1;
  }
  if(maxr>nzt) maxr=nzt;

  // Compute norms
  Cnorm[1]=T1.norm();
  Cnorm[2]=T2.norm();
  Cnorm[3]=T3.norm();
  Cnorm[4]=T4.norm();
  Cnorm[5]=T5.norm();
  Cnorm[6]=T6.norm();
  Cnorm[7]=T7.norm();
  Cnorm[8]=T8.norm();
  Cnorm[9]=T9.norm();
  Cnorm[10]=T10.norm();
  Cnorm[11]=T11.norm();
  Cnorm[12]=T12.norm();
  Cnorm[13]=T13.norm();
  Cnorm[14]=T14.norm();
  Cnorm[15]=T15.norm();
  Cnorm[16]=T16.norm();
#include "evaluate.cc.in"
  Tnorm[1]=T1.norm();
  Tnorm[2]=T2.norm();
  Tnorm[3]=T3.norm();
  Tnorm[4]=T4.norm();
  Tnorm[5]=T5.norm();
  Tnorm[6]=T6.norm();
  Tnorm[7]=T7.norm();
  Tnorm[8]=T8.norm();

  printf("\n%2s %12s %12s %12s\n","n","|C_n|","|T_n|","|T_n|/|C_n|");
  for(int i=1;i<=maxr;i++)
    printf("%2i %e %e %e\n",i,Cnorm[i],Tnorm[i],Tnorm[i]/Cnorm[i]);
  for(int i=maxr+1;i<=16;i++)
    printf("%2i %e\n",i,Cnorm[i]);

  auto tstop(std::chrono::system_clock::now());
  std::chrono::duration<double> elapsed = tstop-tstart;
  printf("Program run in %8.3f seconds\n",elapsed.count());
  fflush(stdout);

  return 0;
}
