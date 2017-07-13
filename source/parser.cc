/*
 Written by Susi Lehtola, 2016
 Copyright (c) 2016, Susi Lehtola

 This program is free software; you can redistribute it and/or
 modify it under the terms of the GNU General Public License
 as published by the Free Software Foundation; either version 2
 of the License, or (at your option) any later version.
*/

#include "parser.hh"
#include <cmath>
#include <cfloat>
#include <cstdlib>
#include <cstdio>
#include <stdexcept>
#include <algorithm>
#include <sstream>
#include <chrono>
#include "wick.hh"

//#define DEBUG

bool operator==(const determinant_t & lhs, const determinant_t & rhs) {
  if(lhs.rank != rhs.rank)
    return false;
  if(lhs.occ.size() != rhs.occ.size())
    throw std::logic_error("Occupied blocks of different size!\n");
  if(lhs.virt.size() != rhs.virt.size())
    throw std::logic_error("Virtual blocks of different size!\n");

  for(size_t i=0;i<lhs.occ.size();i++)
    if(lhs.occ[i] != rhs.occ[i])
      return false;
  for(size_t i=0;i<lhs.virt.size();i++)
    if(lhs.virt[i] != rhs.virt[i])
      return false;

  return true;
}

bool operator<(const determinant_string_t & lhs, const determinant_string_t & rhs) {
  return std::abs(lhs.coeff) > std::abs(rhs.coeff);
}

bool operator<(const determinant_t & lhs, const determinant_t & rhs) {
  if(lhs.rank < rhs.rank)
    return true;
  else if(lhs.rank > rhs.rank)
    return false;

  if(lhs.occ.size() != rhs.occ.size()) {
    std::ostringstream oss;
    oss << "Occupied blocks of different size: " << lhs.occ.size() << " vs " << rhs.occ.size() << "!\n";
    throw std::logic_error(oss.str());
  }
  if(lhs.virt.size() != rhs.virt.size()) {
    std::ostringstream oss;
    oss << "Virtual blocks of different size: " << lhs.virt.size() << " vs " << rhs.virt.size() << "!\n";
    throw std::logic_error(oss.str());
  }

  for(size_t i=0;i<lhs.occ.size();i++)
    if(lhs.occ[i] < rhs.occ[i])
      return true;
    else if(lhs.occ[i] > rhs.occ[i])
      return false;

  for(size_t i=0;i<lhs.virt.size();i++)
    if(lhs.virt[i] < rhs.virt[i])
      return true;
    else if(lhs.virt[i] > rhs.virt[i])
      return false;

  return false;
}

determinant_t operator+(const determinant_t & lhs, const determinant_t & rhs) {
  if(!(lhs==rhs))
    throw std::logic_error("Impossible addition!\n");

  determinant_t ret(lhs);
  ret.coeff+=rhs.coeff;

  return ret;
}

void print(const determinant_t & det, size_t i) {
  printf("Determinant %i, coefficient %e, excitation rank %i:\n",(int) i, det.coeff, det.rank);

  printf("Occupied string: ");
  for(size_t i=0;i<det.occ.size();i++)
    printf(printentry,det.occ[i]);
  printf("\n");
  printf("Virtual  string: ");
  for(size_t i=0;i<det.virt.size();i++)
    printf(printentry,det.virt[i]);
  printf("\n");
}

void check_duplicates(std::vector<determinant_t> dets) {
  // Order the determinant list wrt the strings
  std::sort(dets.begin(),dets.end());

  // Check for duplicate determinants. The list is sorted so we can
  // just check adjacent elements
  for(size_t idet=0;idet+1<dets.size();idet++)
    if(dets[idet] == dets[idet+1]) {
      print(dets[idet],idet+1);
      print(dets[idet+1],idet+2);
      fflush(stdout);

      std::ostringstream oss;
      oss << "Error: determinants " << idet+1 << " and " << idet+2 << " are the same!\n";
      throw std::logic_error(oss.str());
    }

#ifdef DEBUG
  for(size_t jdet=0;jdet<dets.size();jdet++)
    for(size_t idet=0;idet<jdet;idet++)
    if(dets[idet] == dets[jdet]) {
      print(dets[idet],idet+1);
      print(dets[jdet],jdet+1);
      fflush(stdout);

      std::ostringstream oss;
      oss << "Error: determinants " << idet+1 << " and " << jdet+1 << " are the same!\n";
      throw std::logic_error(oss.str());
    }
  printf("Determinants checked for uniqueness.\n");
  fflush(stdout);
#endif
}

std::vector<determinant_t> parse_determinant_strings(const std::vector<determinant_string_t> & strs, size_t Nalpha, size_t Nbeta) {
  std::vector<determinant_t> dets(strs.size());

#ifdef _OPENMP
#pragma omp parallel for
#endif
  for(size_t idet=0;idet<strs.size();idet++) {
    // Alpha and beta occupation strings are
    const std::vector<bool> & aocc(strs[idet].astr);
    const std::vector<bool> & bocc(strs[idet].bstr);
    
    // Count alpha and beta occupations
    size_t oa=0;
    for(size_t io=0;io<aocc.size();io++)
      if(aocc[io])
	oa++;
    size_t va=aocc.size()-oa;
    size_t ob=0;
    for(size_t io=0;io<bocc.size();io++)
      if(bocc[io])
	ob++;
    size_t vb=bocc.size()-ob;

    if(oa!=Nalpha) {
      std::ostringstream oss;
      oss << "Number of alpha electrons " << oa << " in determinant " << idet+1 << " is wrong: it should be " << Nalpha << ".\n";
      //oss << "Alpha string is " << aocc << ".\n";
      throw std::runtime_error(oss.str());
    }
    if(ob!=Nbeta) {
      std::ostringstream oss;
      oss << "Number of beta electrons " << ob << " in determinant " << idet+1 << " is wrong: it should be " << Nbeta << ".\n";
      //oss << "Beta string is " << bocc << ".\n";
      throw std::runtime_error(oss.str());
    }

    // Spin-orbital occupied orbital string (wrt reference) for the
    // configuration
    dets[idet].occ.resize(oa+ob);
    for(size_t io=0;io<oa;io++)
      dets[idet].occ[io]=aocc[io];
    for(size_t io=0;io<ob;io++)
      dets[idet].occ[io+oa]=bocc[io];

    dets[idet].virt.resize(va+vb);
    for(size_t iv=0;iv<va;iv++)
      dets[idet].virt[iv]=aocc[iv+oa];
    for(size_t iv=0;iv<vb;iv++)
      dets[idet].virt[iv+va]=bocc[iv+ob];

    // Determine excitation rank
    size_t nh=0;
    for(size_t io=0;io<dets[idet].occ.size();io++)
      if(!dets[idet].occ[io])
	nh++;
    size_t np=0;
    for(size_t io=0;io<dets[idet].virt.size();io++)
      if(dets[idet].virt[io])
	np++;
    if(nh!=np) throw std::logic_error("Number of particles and holes doesn't match!\n");
    dets[idet].rank=nh;

    // Figure out the sign of the excitation amplitude
    dets[idet].coeff=strs[idet].coeff*string_sign(aocc)*string_sign(bocc);
  }

  return dets;
}

std::vector<determinant_t> parse(const std::string & fname, size_t nmax, double thr) {
  // Read in determinants
  size_t Norb, Nalpha, Nbeta;

  auto tstart(std::chrono::system_clock::now());
  std::vector<determinant_string_t> detstrs(read_determinant_strings(fname,Norb,Nalpha,Nbeta));
  auto tstop(std::chrono::system_clock::now());

  std::chrono::duration<double> elapsed = tstop-tstart;
  printf("%i determinant strings read in %.3f seconds.\n",(int) detstrs.size(),elapsed.count());
  fflush(stdout);
  tstart=tstop;

  // Sort into descending order in weight
  std::sort(detstrs.begin(),detstrs.end());

  // Construct histogram from 1e0 to 1e-16
  std::vector<size_t> hist(17);
  for(size_t i=0;i<detstrs.size();i++) {
    if(std::abs(detstrs[i].coeff)<1e-16) {
      hist[hist.size()-1]++;
    } else {
      double lg(-log10(std::abs(detstrs[i].coeff)));
      int pos(round(lg));
      if(pos<0 || pos>=(int) hist.size())
	throw std::logic_error("Invalid position!\n");
      hist[pos]++;
    }
  }
  // Erase trailing zero elements
  while(hist[hist.size()-1]==0)
    hist.erase(hist.end()-1);
  printf("\nLog10 coefficient histogram\n");
  for(size_t i=0;i<hist.size();i++)
    printf("%2i %i\n",(int) i,(int) hist[i]);
  printf("\n");
  
  // Discard any extra entries
  if(nmax>0) {
    if(detstrs.size()>nmax) {
      printf("%i determinants dropped.\n",(int) (detstrs.size()-nmax));
      detstrs.erase(detstrs.begin()+nmax,detstrs.end());
    }
  } else if(thr>0) {
    for(nmax=0;nmax<detstrs.size();nmax++)
      if(std::abs(detstrs[nmax].coeff)<thr)
	break;
    printf("%i determinants dropped due to small weight, threshold = %e.\n",(int) (detstrs.size()-nmax),thr);
    detstrs.erase(detstrs.begin()+nmax,detstrs.end());
  }

  tstop=std::chrono::system_clock::now();
  elapsed = tstop-tstart;
  printf("Determinant strings sorted and dropped in %.3f seconds.\n",elapsed.count());
  fflush(stdout);
  tstart=tstop;
  
  printf("Running with %i determinants.\n\n",(int) detstrs.size());

  
  std::vector<determinant_t> dets(parse_determinant_strings(detstrs,Nalpha,Nbeta));
  tstop=std::chrono::system_clock::now();
  elapsed = tstop-tstart;
  printf("Determinant strings parsed in %.3f seconds.\n",elapsed.count());
  fflush(stdout);
  
  // Check for duplicates
  check_duplicates(dets);

  // Count distribution of ranks and weights
  std::vector<long unsigned> detrank;
  std::vector<double> detweight;
  for(size_t idet=0;idet<dets.size();idet++) {
    if(dets[idet].rank>=detrank.size()) {
      detrank.resize(dets[idet].rank+1,0);
      detweight.resize(dets[idet].rank+1,0);
    }
    detrank[dets[idet].rank]++;
    detweight[dets[idet].rank]+=dets[idet].coeff*dets[idet].coeff;
  }
  std::vector<long unsigned> cumud(detrank);
  for(size_t i=1;i<cumud.size();i++)
    cumud[i]+=cumud[i-1];
  std::vector<double> cumuw(detweight);
  for(size_t i=1;i<cumuw.size();i++)
    cumuw[i]+=cumuw[i-1];
  printf("%2s %9s %5s %5s %12s %5s\n","n","Ndets","frac","cumu","w","cumu");
  for(size_t r=0;r<detrank.size();r++)
    printf("%2i %9lu %5.1f %5.1f %e %5.1f\n",(int) r, detrank[r], detrank[r]*100.0/dets.size(), cumud[r]*100.0/dets.size(), detweight[r]*100.0, cumuw[r]*100.0);
  printf("\n");

  // Find reference configuration
  std::vector<size_t> refc;
  for(size_t idet=0;idet<dets.size();idet++)
    if(dets[idet].rank==0) {
      refc.push_back(idet);
      print(dets[idet],idet);
      fflush(stdout);
    }
  if(refc.size()!=1) {
    std::ostringstream oss;
    oss << "Invalid reference configuration! Found " << refc.size() << " rank-0 excitations!\n";
    throw std::logic_error(oss.str());
  }

  // Print dominant determinant
  if(refc[0]!=0) {
    printf("Dominant determinant is\n");
    print(dets[0],0);
    fflush(stdout);
  }
  
  // Convert to intermediate normalization
  const double refC(dets[refc[0]].coeff);
  for(size_t idet=0;idet<dets.size();idet++) {
    dets[idet].coeff/=refC;
  }

  return dets;
}

int string_sign(const std::vector<bool> & str) {
  // Count number of electrons
  int nel=0;
  for(size_t i=0;i<str.size();i++)
    if(str[i])
      nel++;

  // Hole and particle indices
  std::vector<int> hs, ps;
  for(int i=0;i<nel;i++)
    if(!str[i])
      hs.push_back(i);
  for(size_t i=nel;i<str.size();i++)
    if(str[i])
      ps.push_back(i);
  if(hs.size() != ps.size())
    throw std::logic_error("Number of holes does not match that of particles!\n");

  // Ket operator string is
  std::vector<int> ketstr;
  for(size_t i=0;i<str.size();i++)
    if(str[i])
      ketstr.push_back(i);

  // while the bra string starts out with the excitation operator
  std::vector<int> brastr(ps);
  for(size_t i=hs.size()-1;i<hs.size();i--)
    brastr.push_back(-hs[i]);
  // and continues by building the HF determinant
  for(int i=0;i<nel;i++)
    brastr.push_back(i);

  // Then, the total operator string is going to be
  std::vector<int> opstr(brastr.size()+ketstr.size());
  for(size_t i=0;i<brastr.size();i++)
    // Hermitian conjugate
    opstr[i]=-brastr[brastr.size()-1-i];
  for(size_t i=0;i<ketstr.size();i++)
    opstr[i+brastr.size()]=ketstr[i];

  // The sign is then given by the Wick contraction
  return wick(opstr);
}
