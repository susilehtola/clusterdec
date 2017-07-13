/*
 Written by Susi Lehtola, 2016
 Copyright (c) 2016, Susi Lehtola

 This program is free software; you can redistribute it and/or
 modify it under the terms of the GNU General Public License
 as published by the Free Software Foundation; either version 2
 of the License, or (at your option) any later version.
*/

#include "wick.hh"
#include <cstdio>
#include <cstdlib>

int wick(const std::vector<int> & opstr) {
  std::vector<size_t> begin, end;
  for(size_t i=0;i<opstr.size();i++)
    for(size_t j=i+1;j<opstr.size();j++)
      if(opstr[i]==-opstr[j]) {
	begin.push_back(i);
	end.push_back(j);
      }

  // Print out contractions
  //for(size_t i=0;i<begin.size();i++)
  //  printf("Contraction %i: begins at %i, ends at %i\n",(int) i,begin[i],end[i]);

  // Loop over products
  int sign=1;
  for(size_t i=0;i<begin.size();i++)
    for(size_t j=0;j<i;j++) {
      if((begin[j]<begin[i] && end[i]<end[j]) || (begin[i]<begin[j] && end[j]<end[i]))
	// Stacked
	continue;
      else if(begin[i]>end[j] || end[j]<begin[i])
	// No touching
	continue;
      else
	// Crossing lines
	sign=-sign;
    }
  
#ifdef DEBUG
  printf("Sign for string <0|");
  for(size_t i=0;i<opstr.size();i++) {
    if(opstr[i]>=0)
      printf(" %i+",opstr[i]);
    else
      printf(" %i-",-opstr[i]);
  }
  printf(" |0> is %i\n",sign);
#endif

  return sign;
}
