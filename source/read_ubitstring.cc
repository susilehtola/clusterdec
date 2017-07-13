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

void convert_to_string(const std::string & num, std::vector<bool> & astr, std::vector<bool> & bstr) {
  // Convert to bit string
  astr.resize(num.size());
  bstr.resize(num.size());
  for(size_t i=0;i<num.size();i++)
    switch(num[i]) {
    case('0'):
      astr[i]=false;
      bstr[i]=false;
      break;
    case('u'):
      astr[i]=true;
      bstr[i]=false;
      break;
    case('d'):
      astr[i]=false;
      bstr[i]=true;
      break;
    case('2'):
      astr[i]=true;
      bstr[i]=true;
      break;
    default:
      std::ostringstream oss;
      oss << "Invalid bit string " << num << "!\n";
      throw std::logic_error(oss.str());
    }
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
  std::vector<std::string> sstr(Ndet);
  for(size_t idet=0;idet<Ndet;idet++) {
    // Get line
    getline(in,line);
    std::istringstream ss(line);

    // Read coefficient and spin string
    ss >> dets[idet].coeff >> sstr[idet];
  }

#ifdef _OPENMP
#pragma omp parallel for
#endif
  for(size_t idet=0;idet<Ndet;idet++)
    // Spin-up and spin-down occupations
    convert_to_string(sstr[idet],dets[idet].astr,dets[idet].bstr);

  return dets;
}
