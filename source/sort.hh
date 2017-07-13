/*
 Written by Susi Lehtola, 2016
 Copyright (c) 2016, Susi Lehtola

 This program is free software; you can redistribute it and/or
 modify it under the terms of the GNU General Public License
 as published by the Free Software Foundation; either version 2
 of the License, or (at your option) any later version.
*/

#ifndef SORT_H
#define SORT_H

#include <cstddef>
#include <vector>
#include "types.hh"

typedef struct {
  /// Original index
  size_t io;
  /// Value
  orbind_t val;
} sort_helper_t;

/// Compare by value
bool operator<(const sort_helper_t & lh, const sort_helper_t & rh);

/// Sort indices into order, also return sign of permutation
std::vector<orbind_t> sort(const std::vector<orbind_t> & in, int & sign);

/// Calculate sign of permutation
int permutation_sign(const std::vector<size_t> & p);

#endif

