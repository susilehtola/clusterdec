/*
 Written by Susi Lehtola, 2016
 Copyright (c) 2016, Susi Lehtola

 This program is free software; you can redistribute it and/or
 modify it under the terms of the GNU General Public License
 as published by the Free Software Foundation; either version 2
 of the License, or (at your option) any later version.
*/

#include "parser.hh"
#include <fstream>
#include <sstream>
#include <stdexcept>

std::vector<bool> convert_to_string(const std::string & num) {
  // Convert to bit string
  std::vector<bool> ret(num.size());
  for(size_t i=0;i<num.size();i++)
    switch(num[i]) {
    case('0'):
      ret[i]=false;
      break;
    case('1'):
      ret[i]=true;
      break;
    default:
      std::ostringstream oss;
      oss << "Invalid bit string " << num << "!\n";
      throw std::logic_error(oss.str());
    }
      
  return ret;
}

std::vector<determinant_string_t> read_determinant_strings(const std::string & fname, size_t & Norb, size_t & Nalpha, size_t & Nbeta) {
  // Open file
  std::ifstream in(fname);
  if(!in.good())
    throw std::runtime_error("Error opening file \"" + fname + "\"!\n");

  std::string line;
  getline(in,line);

  // Stream the line

  // Read number of determinants, orbitals, and alpha and beta electrons
  size_t Ndet;
  {
    std::istringstream ss(line);
    ss >> Ndet >> Norb >> Nalpha >> Nbeta;
  }

  // Read in determinants
  std::vector<determinant_string_t> dets(Ndet);
  std::vector<std::string> astr(Ndet), bstr(Ndet);
  for(size_t idet=0;idet<Ndet;idet++) {
    // Get line
    getline(in,line);
    std::istringstream ss(line);

    // Read coefficient, alpha string and beta string
    ss >> dets[idet].coeff >> astr[idet] >> bstr[idet];
  }

#ifdef _OPENMP
#pragma omp parallel for
#endif
  for(size_t idet=0;idet<Ndet;idet++) {
    // Spin-up and spin-down occupations
    dets[idet].astr=convert_to_string(astr[idet]);
    dets[idet].bstr=convert_to_string(bstr[idet]);
  }

  return dets;
}
