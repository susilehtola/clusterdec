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
// Multiple precision library
#include <gmpxx.h>

std::vector<bool> convert_to_string(mpz_class num, size_t norb) {
  // Convert to bit string
  std::vector<bool> ret;

  mpz_class two("2");
  while(num>0) {
    // num = q*2 + r
    mpz_class b=num%two;
    ret.push_back(b.get_ui());
    num=(num-b)/two;
  }

  // Pad to number of orbitals
  if(ret.size()<norb)
    ret.resize(norb,0);

  // Revert order
  //  std::reverse(ret.begin(), ret.end());

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
  std::vector<mpz_class> astr(Ndet), bstr(Ndet);
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
    if(astr[idet]<0) {
      std::ostringstream oss;
      oss << "Error: got negative value " << astr[idet] << " for alpha string of determinant " << idet+1 << "!\n";
      throw std::runtime_error(oss.str());
    }
    if(bstr[idet]<0) {
      std::ostringstream oss;
      oss << "Error: got negative value " << bstr[idet] << " for beta string of determinant " << idet+1 << "!\n";
      throw std::runtime_error(oss.str());
    }

    // Spin-up and spin-down occupations
    dets[idet].astr=convert_to_string(astr[idet],Norb);
    dets[idet].bstr=convert_to_string(bstr[idet],Norb);

    //printf("Got line \"%s\"; this is a rank-%i determinant.\n",line.c_str(),(int) nh);
  }

  return dets;
}
