/*
  Written by Susi Lehtola, 2016
  Copyright (c) 2016, Susi Lehtola

  This program is free software; you can redistribute it and/or
  modify it under the terms of the GNU General Public License
  as published by the Free Software Foundation; either version 2
  of the License, or (at your option) any later version.
*/

#ifndef GATHER_H
#define GATHER_H

#include <cstdlib>
#include <stdexcept>
#include <vector>
#include "tensor.hh"

inline double gather(const std::vector< std::vector<orbind_t> > & tensors, const std::vector<orbind_t> & oidx, const Tensor1 & T1, const Tensor2 & T2, const Tensor3 & T3, const Tensor4 & T4, const Tensor5 & T5, const Tensor6 & T6, const Tensor7 & T7, const Tensor8 & T8) {
  std::vector<double> tval(tensors.size());

  for(size_t it=0;it<tensors.size();it++) {
    switch(tensors[it].size()) {
    case(2):
      tval[it]=T1(oidx[tensors[it][0]],oidx[tensors[it][1]]);
      break;
    case(4):
      tval[it]=T2(oidx[tensors[it][0]],oidx[tensors[it][1]],oidx[tensors[it][2]],oidx[tensors[it][3]]);
      break;
    case(6):
      tval[it]=T3(oidx[tensors[it][0]],oidx[tensors[it][1]],oidx[tensors[it][2]],oidx[tensors[it][3]],oidx[tensors[it][4]],oidx[tensors[it][5]]);
      break;
    case(8):
      tval[it]=T4(oidx[tensors[it][0]],oidx[tensors[it][1]],oidx[tensors[it][2]],oidx[tensors[it][3]],oidx[tensors[it][4]],oidx[tensors[it][5]],oidx[tensors[it][6]],oidx[tensors[it][7]]);
      break;
    case(10):
      tval[it]=T5(oidx[tensors[it][0]],oidx[tensors[it][1]],oidx[tensors[it][2]],oidx[tensors[it][3]],oidx[tensors[it][4]],oidx[tensors[it][5]],oidx[tensors[it][6]],oidx[tensors[it][7]],oidx[tensors[it][8]],oidx[tensors[it][9]]);
      break;
    case(12):
      tval[it]=T6(oidx[tensors[it][0]],oidx[tensors[it][1]],oidx[tensors[it][2]],oidx[tensors[it][3]],oidx[tensors[it][4]],oidx[tensors[it][5]],oidx[tensors[it][6]],oidx[tensors[it][7]],oidx[tensors[it][8]],oidx[tensors[it][9]],oidx[tensors[it][10]],oidx[tensors[it][11]]);
      break;
    case(14):
      tval[it]=T7(oidx[tensors[it][0]],oidx[tensors[it][1]],oidx[tensors[it][2]],oidx[tensors[it][3]],oidx[tensors[it][4]],oidx[tensors[it][5]],oidx[tensors[it][6]],oidx[tensors[it][7]],oidx[tensors[it][8]],oidx[tensors[it][9]],oidx[tensors[it][10]],oidx[tensors[it][11]],oidx[tensors[it][12]],oidx[tensors[it][13]]);
      break;
    case(16):
      tval[it]=T8(oidx[tensors[it][0]],oidx[tensors[it][1]],oidx[tensors[it][2]],oidx[tensors[it][3]],oidx[tensors[it][4]],oidx[tensors[it][5]],oidx[tensors[it][6]],oidx[tensors[it][7]],oidx[tensors[it][8]],oidx[tensors[it][9]],oidx[tensors[it][10]],oidx[tensors[it][11]],oidx[tensors[it][12]],oidx[tensors[it][13]],oidx[tensors[it][14]],oidx[tensors[it][15]]);
      break;
    default:
      throw std::logic_error("Shouldn't be here!\n");
    }

    // Short-circuit: if element is zero then the whole thing is zero.
    if(tval[it]==0.0)
      return 0.0;
  }

  double dA=1.0;
  for(size_t i=0;i<tval.size();i++)
    dA*=tval[i];
  return dA;
}

#endif
