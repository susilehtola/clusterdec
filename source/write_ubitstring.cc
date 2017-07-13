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
#include <iomanip>

std::string convert_to_string(const std::vector<bool> & astr, const std::vector<bool> & bstr) {
  if(astr.size() != bstr.size())
    throw std::logic_error("Strings are of different size!\n");
  
  // Convert to bit string
  std::ostringstream oss;
  for(size_t i=0;i<astr.size();i++) {
    if(astr[i] && bstr[i])
      oss << '2';
    else if(astr[i] && !bstr[i])
      oss << 'u';
    else if(!astr[i] && bstr[i])
      oss << 'd';
    else if(!astr[i] && !bstr[i])
      oss << '0';
  }
  
  return oss.str();
}

void write_determinant_strings(const std::string & fname, const std::vector<determinant_string_t> & dets, size_t Norb, size_t Nalpha, size_t Nbeta) {
  // Open file
  std::ofstream out(fname);
  if(!out.good())
    throw std::runtime_error("Error opening file \"" + fname + "\"!\n");

  // Write number of determinants, orbitals, and alpha and beta electrons
  out << dets.size() << " " << Norb << " " << Nalpha << " " << Nbeta << std::endl;

  // Write out determinants
  for(size_t idet=0;idet<dets.size();idet++)
    out << std::scientific << std::setprecision(15) << std::setw(22) << dets[idet].coeff << " " << convert_to_string(dets[idet].astr,dets[idet].bstr) << std::endl;
}
