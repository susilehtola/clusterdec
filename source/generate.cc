/*
 Written by Susi Lehtola, 2016
 Copyright (c) 2016, Susi Lehtola

 This program is free software; you can redistribute it and/or
 modify it under the terms of the GNU General Public License
 as published by the Free Software Foundation; either version 2
 of the License, or (at your option) any later version.
*/

#include <cstdio>
#include <cstdlib>
#include <vector>
#include <cmath>
#include <string>
#include <sstream>
#include <algorithm>
#include <chrono>
#include <stdexcept>
#include "wick.hh"

/*
   To figure out which products of excitation tensors can couple up to
   a rank r, we generate an index vector that contains the number of
   T_n tensors.

   In C_n, we can have contributions from T_1 to T_{n-1} (T_n is
   included in the definition!). Each tensor can appear from 0 to n
   times. Thus, we have a n-1 element vector, with indices that have
   n+1 allowed values.

   This function generates the index vector from the product index.
*/
void index(size_t iprod, size_t r, std::vector<size_t> & idx) {
  idx.resize(r-1);

  // Fill out the index
  size_t i=0, quotient(iprod);
  while(quotient>0) {
    // Compute index
    idx[i]=(quotient%(r+1));
    // Adjust quotient
    quotient/=r+1;
    // Increment index
    i++;
  }
}

// Generates list of products of T_n that yield C_r
std::vector< std::vector<size_t> > target_T(size_t r) {
  std::vector< std::vector<size_t> > ret;

  // Loop over all possible products
  for(size_t iprod=0;iprod<std::pow(r+1,r-1);iprod++) {
    std::vector<size_t> idx;
    index(iprod,r,idx);

    // Calculate total rank
    size_t rank=0;
    for(size_t i=0;i<idx.size();i++)
      rank+=(i+1)*idx[i];

    // This is one we want
    if(rank==r)
      ret.push_back(idx);
  }

  /*
  // Debug for checking consistency
  ret.resize(1);
  ret[0].assign(r,0);
  ret[0][r-1]=1;
  */

  return ret;
}

// Amplitude tensor
typedef struct {
  /// Occupied indices
  std::vector<unsigned char> occ;
  /// Virtual indices
  std::vector<unsigned char> virt;
} tensor_t;

/* Sorts the occupied and virtual strings into numerical order. The
   sign doesn't matter, since it's only computed for the canonical
   ordering. */
void sort(tensor_t & t) {
  std::sort(t.occ.begin(),t.occ.end());
  std::sort(t.virt.begin(),t.virt.end());
}

// Comparison for sorts
bool operator<(const tensor_t & lh, const tensor_t & rh) {
  // First, sort by rank
  if(lh.occ.size() < rh.occ.size())
    return true;
  else if(lh.occ.size() > rh.occ.size())
    return false;

  // then by occupied indices
  for(size_t i=0;i<lh.occ.size();i++)
    if(lh.occ[i] < rh.occ[i])
      return true;
    else if(lh.occ[i] > rh.occ[i])
      return false;

  // and then by virtual indices
  for(size_t i=0;i<lh.virt.size();i++)
    if(lh.virt[i] < rh.virt[i])
      return true;
    else if(lh.virt[i] > rh.virt[i])
      return false;

  // Default return
  return false;
}

// Comparison for storage of nonredudant strings
bool operator==(const tensor_t & lh, const tensor_t & rh) {
  if(lh.occ.size() != rh.occ.size())
    return false;
  for(size_t i=0;i<lh.occ.size();i++)
    if(lh.occ[i] != rh.occ[i])
      return false;
  for(size_t i=0;i<lh.virt.size();i++)
    if(lh.virt[i] != rh.virt[i])
      return false;
  return true;
}

// Comparison operator expressed in previously implemented operators
bool operator>(const tensor_t & lh, const tensor_t & rh) {
  return (!(lh<rh) && !(lh==rh));
}

// Comparison operator expressed in previously implemented operators
bool operator!=(const tensor_t & lh, const tensor_t & rh) {
  return !(lh==rh);
}

/// Tensor contraction
typedef struct {
  /// Contraction
  std::vector<tensor_t> c;
  /// Sign
  int sign;
} tensorc_t;

/// Get rank of contraction
int rank(const tensorc_t & t) {
  int r=0;
  for(size_t i=0;i<t.c.size();i++)
    r+=t.c[i].occ.size();
  return r;
}

/**
   Sorts the contraction. We don't care about the sign here, since the
   point is to convert the generated indices into the canonical
   ordering, for which the sign is then computed using Wick
   contractions.
*/
void sort(tensorc_t & t, bool tsort=true) {
  for(size_t i=0;i<t.c.size();i++)
    sort(t.c[i]);
  if(tsort)
    std::stable_sort(t.c.begin(),t.c.end());
}

// Comparison for sorts
bool operator<(const tensorc_t & lh, const tensorc_t & rh) {
  if(lh.c.size() < rh.c.size())
    return true;
  else if(lh.c.size() > rh.c.size())
    return false;

  // Compare individual terms
  for(size_t i=0;i<lh.c.size();i++)
    if(lh.c[i] < rh.c[i])
      return true;
    else if(lh.c[i] > rh.c[i])
      return false;

  return false;
}

// Comparison for nonredundant storage
bool operator==(const tensorc_t & lh, const tensorc_t & rh) {
  if(lh.c.size() != rh.c.size())
    return false;
  for(size_t i=0;i<lh.c.size();i++)
    if(lh.c[i] != rh.c[i])
      return false;
  return true;
}

bool operator!=(const tensorc_t & lh, const tensorc_t & rh) {
  return !(lh==rh);
}

void print(const tensor_t & t) {
  if(t.occ.size())
    printf("t(");
  for(size_t i=0;i<t.occ.size();i++)
    printf(" i%i",t.occ[i]);

  for(size_t i=0;i<t.virt.size();i++)
    printf(" a%i",t.virt[i]);
  if(t.virt.size())
    printf(" )");
}

void print(const tensorc_t & t) {
  if(t.sign==1)
    printf("+");
  else printf("-");

  for(size_t i=0;i<t.c.size();i++) {
    if(i)
      printf(" ");
    print(t.c[i]);
  }
}

void calc_sign(tensorc_t & t) {
  // Get the sign from the Wick contraction.

  int r(rank(t));

  // Operator string
  std::vector<int> opstr(4*r);

  // The hole creation operators come first in order
  for(int i=0;i<r;i++)
    opstr[i]=i+1;
  // after which follow the particle annihilation indices in reverse order
  for(int i=0;i<r;i++)
    opstr[i+r]=-(2*r-i);

  // Now, we put in the string indices
  {
    size_t idx=2*r;
    for(size_t ic=0;ic<t.c.size();ic++) {
      // First come the creation indices
      for(size_t iv=0;iv<t.c[ic].virt.size();iv++)
	opstr[idx++]=t.c[ic].virt[iv]+1+r;
      // and then the hole annihilation indices in opposite order
      for(size_t ih=t.c[ic].occ.size()-1;ih<t.c[ic].occ.size();ih--)
	opstr[idx++]=-t.c[ic].occ[ih]-1;
    }
  }

  // Full list of contractions is then
  t.sign=wick(opstr);
}

// Pair tensors
std::vector<tensorc_t> form_pairings(size_t r) {
  // Returned array
  std::vector<tensorc_t> ret;

  // Get list of tensors that pair to rank r
  std::vector< std::vector<size_t> > prod(target_T(r));

  // Loop over the combinations
  for(size_t icomb=0;icomb<prod.size();icomb++) {
    auto tstart(std::chrono::system_clock::now());

    {
      std::ostringstream oss;
      for(size_t i=0;i<prod[icomb].size();i++)
	if(prod[icomb][i]) {
	  oss << " T" << i+1;
	  if(prod[icomb][i]>1)
	    oss << "^" << prod[icomb][i];
	}
      printf("%40s: ... ",oss.str().c_str());
      fflush(stdout);
    }

    // Generate the tensors in the canonical order
    tensorc_t orig;
    {
      // Total rank so far
      int rr=0;

      // Loop over tensor rank
      for(size_t it=0;it<prod[icomb].size();it++) {
	// Loop over appearances of the tensor
	for(size_t ia=0;ia<prod[icomb][it];ia++) {
	  tensor_t add;

	  // Add occupied indices
	  for(size_t ii=0;ii<=it;ii++) // it is off by one in rank
	    add.occ.push_back(ii+rr);
	  // Add virtual indices
	  for(size_t ii=0;ii<=it;ii++)
	    add.virt.push_back(ii+rr);
	  // Increment counter
	  rr+=it+1;

	  // Add tensor
	  orig.c.push_back(add);
	}
      }
    }
    // Dummy sign
    orig.sign=0;
    // Sort
    sort(orig);

    // Indices to permute
    std::vector<int> idx(r);
    for(size_t i=0;i<r;i++)
      idx[i]=i;

    // Because we can permute the indices within any given tensor,
    // the generation takes a long time if we loop over all possible
    // permutations. E.g. the T1*T7 term in the octuples has (8!)^2
    // permutations out of which (7!)^2 correspond to permutations
    // inside T7.
    std::vector< std::vector<int> > permgen;
    {
      std::vector<tensorc_t> sperms;
      {
	std::vector<int> iidx(idx);
	do {
	  // New permuted contraction is
	  tensorc_t pt(orig);
	  // Apply index permutation to occupied indices
	  for(size_t ic=0;ic<pt.c.size();ic++)
	    for(size_t io=0;io<pt.c[ic].occ.size();io++)
	      pt.c[ic].occ[io]=iidx[orig.c[ic].occ[io]];
	  for(size_t ic=0;ic<pt.c.size();ic++)
	    for(size_t io=0;io<pt.c[ic].virt.size();io++)
	      pt.c[ic].virt[io]=iidx[orig.c[ic].virt[io]];
	  // Sort indices, not tensors
	  sort(pt,false);

	  // Check if an equivalent contraction is on the list
	  std::vector<tensorc_t>::iterator low(std::lower_bound(sperms.begin(),sperms.end(),pt));
	  if(low == sperms.end() || *low != pt) {
	    sperms.insert(low,pt);
	    permgen.push_back(iidx);
	  }
	} while(std::next_permutation(iidx.begin(),iidx.end()));
      }
    }

    // Now we can loop over the reduced set of permutations with ease,
    // in parallel!
    std::vector<tensorc_t> perms;

#ifdef _OPENMP
#pragma omp parallel
    {
      // Private copy of permutations
      std::vector<tensorc_t> pperms;

      size_t N(permgen.size());
      size_t N2(N*N);

      /*
#pragma omp for collapse(2)
      for(size_t ip=0;ip<permgen.size();ip++)
	for(size_t jp=0;jp<permgen.size();jp++) {
      */
#pragma omp for schedule(dynamic)
      for(size_t ijp=0;ijp<N2;ijp++) {
	size_t ip=ijp%N;
	size_t jp=ijp/N;
	  // New permuted contraction is
	  tensorc_t pt(orig);
	  // Apply index permutation to occupied indices
	  for(size_t ic=0;ic<pt.c.size();ic++)
	    for(size_t io=0;io<pt.c[ic].occ.size();io++)
	      pt.c[ic].occ[io]=permgen[ip][pt.c[ic].occ[io]];
	  for(size_t ic=0;ic<pt.c.size();ic++)
	    for(size_t io=0;io<pt.c[ic].virt.size();io++)
	      pt.c[ic].virt[io]=permgen[jp][pt.c[ic].virt[io]];
	  // Sort
	  sort(pt);

	  // Check if an equivalent contraction is on the list
	  std::vector<tensorc_t>::iterator low(std::lower_bound(pperms.begin(),pperms.end(),pt));
	  if(low == pperms.end() || *low != pt) {
	    // Calculate sign
	    calc_sign(pt);
	    // and add to list
	    pperms.insert(low,pt);
	  }
	}

#pragma omp critical
      {
	// Collect results: loop over private copies
	for(size_t i=0;i<pperms.size();i++) {
	  tensorc_t pt(pperms[i]);

	  // and add them to the full list
	  std::vector<tensorc_t>::iterator low(std::lower_bound(perms.begin(),perms.end(),pt));
	  if(low == perms.end() || *low != pt) {
	    // and add to list if they don't exist
	    perms.insert(low,pt);
	  }
	}
      }
    }

#else
    for(size_t ip=0;ip<permgen.size();ip++)
      for(size_t jp=0;jp<permgen.size();jp++) {
	// New permuted contraction is
	tensorc_t pt(orig);
	// Apply index permutation to occupied indices
	for(size_t ic=0;ic<pt.c.size();ic++)
	  for(size_t io=0;io<pt.c[ic].occ.size();io++)
	    pt.c[ic].occ[io]=permgen[ip][pt.c[ic].occ[io]];
	for(size_t ic=0;ic<pt.c.size();ic++)
	  for(size_t io=0;io<pt.c[ic].virt.size();io++)
	    pt.c[ic].virt[io]=permgen[jp][pt.c[ic].virt[io]];
	// Sort
	sort(pt);

	// Check if an equivalent contraction is on the list
	std::vector<tensorc_t>::iterator low(std::lower_bound(perms.begin(),perms.end(),pt));
	if(low == perms.end() || *low != pt) {
	  // Calculate sign
	  calc_sign(pt);
	  // and add to list
	  perms.insert(low,pt);
	}
      }
#endif

    ret.insert(ret.end(),perms.begin(),perms.end());

    auto tstop(std::chrono::system_clock::now());
    std::chrono::duration<double> elapsed = tstop-tstart;
    fprintf(stderr,"%8lu terms (%9.3f s)\n",(long unsigned) perms.size(),elapsed.count());
  }

  return ret;
}

void print_code(const std::vector<tensorc_t> & list, FILE * out) {
  int r(rank(list[0]));

  fprintf(out,"  /* <%i|exp(T)|%i> */\n#ifdef _OPENMP\n#pragma omp parallel for\n#endif\n",r,r);
  fprintf(out,"  for(size_t idet=0;idet<T%i.size();idet++) {\n",r);
  fprintf(out,"    Tensor%i_entry_t etr(T%i.get(idet));\n",r,r);

  fprintf(out,"    int");
  for(int i=0;i<r;i++) {
    if(i)
      fprintf(out,",");
    fprintf(out," o%i",i);
  }
  fprintf(out,";\n");

  fprintf(out,"    int");
  for(int i=0;i<r;i++) {
    if(i)
      fprintf(out,",");
    fprintf(out," v%i",i);
  }
  fprintf(out,";\n");

  fprintf(out,"    double A;\n");
  fprintf(out,"    unroll_entry%i(etr",r);
  for(int i=0;i<r;i++)
    fprintf(out,", o%i",i);
  for(int i=0;i<r;i++)
    fprintf(out,", v%i",i);
  fprintf(out,", A);\n");

  for(size_t i=0;i<list.size();i++) {
    // Take sign into account here
    fprintf(out,"    A += %+1i",-list[i].sign);

    for(size_t ic=0;ic<list[i].c.size();ic++) {
      fprintf(out,"*T%i(",(int) list[i].c[ic].occ.size());
      for(size_t io=0;io<list[i].c[ic].occ.size();io++) {
	if(io)
	  fprintf(out,",");
	fprintf(out,"o%i",list[i].c[ic].occ[io]);
      }
      for(size_t iv=0;iv<list[i].c[ic].virt.size();iv++) {
	fprintf(out,",v%i",list[i].c[ic].virt[iv]);
      }
      fprintf(out,")");
    }
    fprintf(out,";\n");
  }

  fprintf(out,"    etr.el = A;\n");
  fprintf(out,"    T%i.set(idet,etr);\n",r);
  fprintf(out,"  }\n");
  fprintf(out,"  printf(\"%i %%e\\n\",T%i.norm();\n",r,r);
  fflush(out);
}

void print_recipe(const std::vector<tensorc_t> & list) {
  if(!list.size())
    throw std::logic_error("Empty recipe list!\n");

  int r(rank(list[0]));
  char fname[80];
  sprintf(fname,"recipes_%i.dat",r);
  FILE *out=fopen(fname,"w");

  for(size_t i=0;i<list.size();i++) {
    // Take sign into account here
    fprintf(out,"%+i %i",-list[i].sign,(int) list[i].c.size());

    for(size_t ic=0;ic<list[i].c.size();ic++) {
      fprintf(out," %i",(int) list[i].c[ic].occ.size());
      for(size_t io=0;io<list[i].c[ic].occ.size();io++)
	fprintf(out," %i",list[i].c[ic].occ[io]);
      // Shift virtual indices by rank
      for(size_t iv=0;iv<list[i].c[ic].virt.size();iv++)
	fprintf(out," %i",list[i].c[ic].virt[iv]+r);
    }

    fprintf(out,"\n");
  }
  fclose(out);
}

int main(int argc, char **argv) {
  if(argc>2) {
    printf("Usage: %s (rank)\n",argv[0]);
    return 1;
  }

  if(argc==1) {
    FILE *out=fopen("main.cc","w");
    fprintf(out,"/* THIS FILE IS AUTOGENERATED, DON'T TOUCH IT! */\n\n");
    fprintf(out,"#include \"gather.hh\"\n");
    fprintf(out,"#include \"parser.hh\"\n");
    fprintf(out,"#include \"tensor.hh\"\n");
    fprintf(out,"#include \"types.hh\"\n");
    fprintf(out,"#include <vector>\n");
    fprintf(out,"#include <climits>\n");
    fprintf(out,"#include <cmath>\n");
    fprintf(out,"#include <chrono>\n");
    fprintf(out,"#include <sstream>\n");
    fprintf(out,"#include <cstdio>\n\n");

    fprintf(out,"/*\n Written by Susi Lehtola, 2016\n Copyright (c) 2016, Susi Lehtola\n\n This program is free software; you can redistribute it and/or\n modify it under the terms of the GNU General Public License\n as published by the Free Software Foundation; either version 2\n of the License, or (at your option) any later version.\n*/\n");
    fprintf(out,"/* This code is autogenerated by generate.cc. Don't edit it by hand. */\n");

    // Go up to rank
    int maxr=9;
    // Tensors maximum
    int tmax=16;

    fprintf(out,"int main(int argc, char **argv) {\n");
    fprintf(out,"  printf(\"CLUSTERDEC, written by Susi Lehtola (2016)\\n\\n\");\n");
    fprintf(out,"  if(argc>4) {\n");
    fprintf(out,"    printf(\"Usage: %%s (file) (maxr) (ndets)\\n\",argv[0]);\n");
    fprintf(out,"    return 1;\n");
    fprintf(out,"  }\n");
    fprintf(out,"  auto tstart(std::chrono::system_clock::now());\n");
    fprintf(out,"  std::string fname = (argc>1) ? argv[1] : \"wfs.dat\";\n");
    fprintf(out,"  int maxr = (argc>2) ? atoi(argv[2]) : %i;\n",maxr);
    fprintf(out,"  int ndets = 0;\n");
    fprintf(out,"  double thr = 1e-8;\n");
    fprintf(out,"  if(argc>3) {\n");
    fprintf(out,"    if(atof(argv[3])>1) {\n");
    fprintf(out,"      ndets=atoi(argv[3]);\n");
    fprintf(out,"    } else {\n");
    fprintf(out,"      thr=atof(argv[3]);\n");
    fprintf(out,"    }\n");
    fprintf(out,"  }\n");
    fprintf(out,"  std::vector<determinant_t> list(parse(fname,ndets,thr));\n");
    fprintf(out,"  std::vector<double> Cnorm(%i,0.0), Tnorm(%i,0.0);\n",tmax+1,maxr+1);

    fprintf(out,"\n  // Form excitation tensors\n");
    for(int i=1;i<=tmax;i++) {
      fprintf(out,"  Tensor%i T%i(list);\n",i,i);
    }

    fprintf(out,"\n  // Check maximum excitation rank\n");
    fprintf(out,"  int nzt=0;\n");
    for(int i=maxr;i>=1;i--) {
      if(i==maxr)
	fprintf(out,"  if(T%i.size()) {\n",i);
      else
	fprintf(out,"  } else if(T%i.size()) {\n",i);
      fprintf(out,"      nzt=%i;\n",i);
    }
    fprintf(out,"  }\n");
    fprintf(out,"  if(maxr>nzt) maxr=nzt;\n");

    fprintf(out,"\n  // Compute norms\n");
    for(int i=1;i<=tmax;i++)
      fprintf(out,"  Cnorm[%i]=T%i.norm();\n",i,i);

    fprintf(out,"#include \"evaluate.cc.in\"\n");

    for(int r=1;r<=maxr;r++)
      fprintf(out,"  Tnorm[%i]=T%i.norm();\n",r,r);

    fprintf(out,"\n  printf(\"\\n%%2s %%12s %%12s %%12s\\n\",\"n\",\"|C_n|\",\"|T_n|\",\"|T_n|/|C_n|\");\n");
    fprintf(out,"  for(int i=1;i<=maxr;i++)\n");
    fprintf(out,"    printf(\"%%2i %%e %%e %%e\\n\",i,Cnorm[i],Tnorm[i],Tnorm[i]/Cnorm[i]);\n");
    fprintf(out,"  for(int i=maxr+1;i<=%i;i++)\n",tmax);
    fprintf(out,"    printf(\"%%2i %%e\\n\",i,Cnorm[i]);\n");

    fprintf(out,"\n  auto tstop(std::chrono::system_clock::now());\n");
    fprintf(out,"  std::chrono::duration<double> elapsed = tstop-tstart;\n");
    fprintf(out,"  printf(\"Program run in %%8.3f seconds\\n\",elapsed.count());\n");
    fprintf(out,"  fflush(stdout);\n\n");
    fprintf(out,"  return 0;\n}\n");
    fclose(out);

  } else {
    int r(atoi(argv[1]));

    // Form the pairings
    auto tstart(std::chrono::system_clock::now());
    std::vector<tensorc_t> list(form_pairings(r));
    auto tstop(std::chrono::system_clock::now());

    std::chrono::duration<double> elapsed = tstop-tstart;
    fprintf(stderr,"rank-%i: %8lu terms (%9.3f s)\n",r,(long unsigned) list.size(),elapsed.count());
    fflush(stderr);

    print_recipe(list);
  }

  return 0;
}
