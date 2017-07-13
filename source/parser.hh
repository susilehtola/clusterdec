/*
 Written by Susi Lehtola, 2016
 Copyright (c) 2016, Susi Lehtola

 This program is free software; you can redistribute it and/or
 modify it under the terms of the GNU General Public License
 as published by the Free Software Foundation; either version 2
 of the License, or (at your option) any later version.
*/

#ifndef PARSER_H
#define PARSER_H

#include <vector>
#include <string>
#include <cinttypes>
#include <cstdint>
#include "types.hh"

/// Number type for orbital numberings
#define entry unsigned long int
/// Print type
#define printentry "%lu"

typedef struct {
  // Alpha string
  std::vector<bool> astr;
  // Beta string
  std::vector<bool> bstr;
  // Coefficient
  double coeff;
} determinant_string_t;

/// Compare weights
bool operator<(const determinant_string_t & lhs, const determinant_string_t & rhs);

typedef struct {
  // Occupied string
  std::vector<entry> occ;
  // Virtual string
  std::vector<entry> virt;
  // Coefficient
  double coeff;
  // Rank of excitation
  unsigned int rank;
} determinant_t;

// Are determinants the same?
bool operator==(const determinant_t & lhs, const determinant_t & rhs);
// Compare strings
bool operator<(const determinant_t & lhs, const determinant_t & rhs);
// Addition
determinant_t operator+(const determinant_t & lhs, const determinant_t & rhs);

void print(const determinant_t & det, size_t i);

// Read in determinants (several implementations)
std::vector<determinant_string_t> read_determinant_strings(const std::string & fname, size_t & Norb, size_t & Nalpha, size_t & Nbeta);
// Write out determinant strings (several implementations, for file conversion)
void write_determinant_strings(const std::string & fname, const std::vector<determinant_string_t> & dets, size_t Norb, size_t Nalpha, size_t Nbeta);

// Parse determinants
std::vector<determinant_t> parse(const std::string & fname, size_t ndets, double thr);

// Calculate the sign for the string
int string_sign(const std::vector<bool> & occ);

/// Recipe type
typedef struct {
  /// Sign of term
  int sign;
  /// Tensors
  std::vector< std::vector<orbind_t> > tensors;
} recipe_t;

#endif
