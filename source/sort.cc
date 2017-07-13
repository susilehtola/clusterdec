/*
 Written by Susi Lehtola, 2016
 Copyright (c) 2016, Susi Lehtola

 This program is free software; you can redistribute it and/or
 modify it under the terms of the GNU General Public License
 as published by the Free Software Foundation; either version 2
 of the License, or (at your option) any later version.
*/

#include "sort.hh"
#include <algorithm>

bool operator<(const sort_helper_t & lh, const sort_helper_t & rh) {
  return lh.val	< rh.val;
}

std::vector<orbind_t> sort(const std::vector<orbind_t> & in, int & sign) {
  // Index array
  std::vector<sort_helper_t> hlp(in.size());
  for(size_t i=0;i<in.size();i++) {
    hlp[i].io=i;
    hlp[i].val=in[i];
  }
  std::sort(hlp.begin(),hlp.end());

  // Permutation that makes the indices in order is
  std::vector<size_t> perm(in.size());
  for(size_t i=0;i<in.size();i++)
    perm[i]=hlp[i].io;
  // which yields the sign
  sign=permutation_sign(perm);

  // Indices in the correct order are
  std::vector<orbind_t> ret(in.size());
  for(size_t i=0;i<in.size();i++) {
    ret[i]=hlp[i].val;
  }

  return ret;
}

int permutation_sign(const std::vector<size_t> & p) {
  std::vector<bool> visited(p.size(),false);
  int sign=1;

  for(size_t k=0;k<p.size();k++) {
    if(!visited[k]) {
      size_t next=k;
      size_t length=0;
      while(!visited[next]) {
        // New cycle, count its length
        length++;
        visited[next]=true;
        next=p[next];
      }
      // If length is even, flip sign
      if(length%2==0)
        sign=-sign;
    }
  }

  return sign;
}
