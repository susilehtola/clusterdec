/*
 Written by Susi Lehtola, 2016
 Copyright (c) 2016, Susi Lehtola

 This program is free software; you can redistribute it and/or
 modify it under the terms of the GNU General Public License
 as published by the Free Software Foundation; either version 2
 of the License, or (at your option) any later version.
*/

/* This code is autogenerated by gentensor.py. Don't edit it by hand. */

#ifndef TENSOR_HH
#define TENSOR_HH
#include <vector>
#include <cstddef>
#include "types.hh"
#include "parser.hh"

typedef struct {
  orbind_t o0;
  orbind_t v0;
  double el;
} Tensor1_entry_t;
bool operator<(const Tensor1_entry_t & lh, const Tensor1_entry_t & rh);
bool operator==(const Tensor1_entry_t & lh, const Tensor1_entry_t & rh);
Tensor1_entry_t make_entry1(orbind_t o0, orbind_t v0);
void unroll_entry1(const Tensor1_entry_t & etr, orbind_t & o0, orbind_t & v0, double & val);
void unroll_entry1(const Tensor1_entry_t & etr, std::vector<orbind_t> & idx, double & val);
class Tensor1 {
  std::vector<Tensor1_entry_t> list;
public:
  Tensor1();
  Tensor1(const std::vector<determinant_t> & list);
  ~Tensor1();
  void add(const Tensor1_entry_t & key);
  size_t size() const;
  double operator()(orbind_t o0, orbind_t v0) const;
  Tensor1_entry_t get(size_t idx) const;
  void set(size_t idx, Tensor1_entry_t val);
  void set(size_t idx, double val);
  double norm() const;
  void print(double thr=1e-10) const;
};
typedef struct {
  orbind_t o0, o1;
  orbind_t v0, v1;
  double el;
} Tensor2_entry_t;
bool operator<(const Tensor2_entry_t & lh, const Tensor2_entry_t & rh);
bool operator==(const Tensor2_entry_t & lh, const Tensor2_entry_t & rh);
Tensor2_entry_t make_entry2(orbind_t o0, orbind_t o1, orbind_t v0, orbind_t v1);
void unroll_entry2(const Tensor2_entry_t & etr, orbind_t & o0, orbind_t & o1, orbind_t & v0, orbind_t & v1, double & val);
void unroll_entry2(const Tensor2_entry_t & etr, std::vector<orbind_t> & idx, double & val);
class Tensor2 {
  std::vector<Tensor2_entry_t> list;
public:
  Tensor2();
  Tensor2(const std::vector<determinant_t> & list);
  ~Tensor2();
  void add(const Tensor2_entry_t & key);
  size_t size() const;
  double operator()(orbind_t o0, orbind_t o1, orbind_t v0, orbind_t v1) const;
  Tensor2_entry_t get(size_t idx) const;
  void set(size_t idx, Tensor2_entry_t val);
  void set(size_t idx, double val);
  double norm() const;
  void print(double thr=1e-10) const;
};
typedef struct {
  orbind_t o0, o1, o2;
  orbind_t v0, v1, v2;
  double el;
} Tensor3_entry_t;
bool operator<(const Tensor3_entry_t & lh, const Tensor3_entry_t & rh);
bool operator==(const Tensor3_entry_t & lh, const Tensor3_entry_t & rh);
Tensor3_entry_t make_entry3(orbind_t o0, orbind_t o1, orbind_t o2, orbind_t v0, orbind_t v1, orbind_t v2);
void unroll_entry3(const Tensor3_entry_t & etr, orbind_t & o0, orbind_t & o1, orbind_t & o2, orbind_t & v0, orbind_t & v1, orbind_t & v2, double & val);
void unroll_entry3(const Tensor3_entry_t & etr, std::vector<orbind_t> & idx, double & val);
class Tensor3 {
  std::vector<Tensor3_entry_t> list;
public:
  Tensor3();
  Tensor3(const std::vector<determinant_t> & list);
  ~Tensor3();
  void add(const Tensor3_entry_t & key);
  size_t size() const;
  double operator()(orbind_t o0, orbind_t o1, orbind_t o2, orbind_t v0, orbind_t v1, orbind_t v2) const;
  Tensor3_entry_t get(size_t idx) const;
  void set(size_t idx, Tensor3_entry_t val);
  void set(size_t idx, double val);
  double norm() const;
  void print(double thr=1e-10) const;
};
typedef struct {
  orbind_t o0, o1, o2, o3;
  orbind_t v0, v1, v2, v3;
  double el;
} Tensor4_entry_t;
bool operator<(const Tensor4_entry_t & lh, const Tensor4_entry_t & rh);
bool operator==(const Tensor4_entry_t & lh, const Tensor4_entry_t & rh);
Tensor4_entry_t make_entry4(orbind_t o0, orbind_t o1, orbind_t o2, orbind_t o3, orbind_t v0, orbind_t v1, orbind_t v2, orbind_t v3);
void unroll_entry4(const Tensor4_entry_t & etr, orbind_t & o0, orbind_t & o1, orbind_t & o2, orbind_t & o3, orbind_t & v0, orbind_t & v1, orbind_t & v2, orbind_t & v3, double & val);
void unroll_entry4(const Tensor4_entry_t & etr, std::vector<orbind_t> & idx, double & val);
class Tensor4 {
  std::vector<Tensor4_entry_t> list;
public:
  Tensor4();
  Tensor4(const std::vector<determinant_t> & list);
  ~Tensor4();
  void add(const Tensor4_entry_t & key);
  size_t size() const;
  double operator()(orbind_t o0, orbind_t o1, orbind_t o2, orbind_t o3, orbind_t v0, orbind_t v1, orbind_t v2, orbind_t v3) const;
  Tensor4_entry_t get(size_t idx) const;
  void set(size_t idx, Tensor4_entry_t val);
  void set(size_t idx, double val);
  double norm() const;
  void print(double thr=1e-10) const;
};
typedef struct {
  orbind_t o0, o1, o2, o3, o4;
  orbind_t v0, v1, v2, v3, v4;
  double el;
} Tensor5_entry_t;
bool operator<(const Tensor5_entry_t & lh, const Tensor5_entry_t & rh);
bool operator==(const Tensor5_entry_t & lh, const Tensor5_entry_t & rh);
Tensor5_entry_t make_entry5(orbind_t o0, orbind_t o1, orbind_t o2, orbind_t o3, orbind_t o4, orbind_t v0, orbind_t v1, orbind_t v2, orbind_t v3, orbind_t v4);
void unroll_entry5(const Tensor5_entry_t & etr, orbind_t & o0, orbind_t & o1, orbind_t & o2, orbind_t & o3, orbind_t & o4, orbind_t & v0, orbind_t & v1, orbind_t & v2, orbind_t & v3, orbind_t & v4, double & val);
void unroll_entry5(const Tensor5_entry_t & etr, std::vector<orbind_t> & idx, double & val);
class Tensor5 {
  std::vector<Tensor5_entry_t> list;
public:
  Tensor5();
  Tensor5(const std::vector<determinant_t> & list);
  ~Tensor5();
  void add(const Tensor5_entry_t & key);
  size_t size() const;
  double operator()(orbind_t o0, orbind_t o1, orbind_t o2, orbind_t o3, orbind_t o4, orbind_t v0, orbind_t v1, orbind_t v2, orbind_t v3, orbind_t v4) const;
  Tensor5_entry_t get(size_t idx) const;
  void set(size_t idx, Tensor5_entry_t val);
  void set(size_t idx, double val);
  double norm() const;
  void print(double thr=1e-10) const;
};
typedef struct {
  orbind_t o0, o1, o2, o3, o4, o5;
  orbind_t v0, v1, v2, v3, v4, v5;
  double el;
} Tensor6_entry_t;
bool operator<(const Tensor6_entry_t & lh, const Tensor6_entry_t & rh);
bool operator==(const Tensor6_entry_t & lh, const Tensor6_entry_t & rh);
Tensor6_entry_t make_entry6(orbind_t o0, orbind_t o1, orbind_t o2, orbind_t o3, orbind_t o4, orbind_t o5, orbind_t v0, orbind_t v1, orbind_t v2, orbind_t v3, orbind_t v4, orbind_t v5);
void unroll_entry6(const Tensor6_entry_t & etr, orbind_t & o0, orbind_t & o1, orbind_t & o2, orbind_t & o3, orbind_t & o4, orbind_t & o5, orbind_t & v0, orbind_t & v1, orbind_t & v2, orbind_t & v3, orbind_t & v4, orbind_t & v5, double & val);
void unroll_entry6(const Tensor6_entry_t & etr, std::vector<orbind_t> & idx, double & val);
class Tensor6 {
  std::vector<Tensor6_entry_t> list;
public:
  Tensor6();
  Tensor6(const std::vector<determinant_t> & list);
  ~Tensor6();
  void add(const Tensor6_entry_t & key);
  size_t size() const;
  double operator()(orbind_t o0, orbind_t o1, orbind_t o2, orbind_t o3, orbind_t o4, orbind_t o5, orbind_t v0, orbind_t v1, orbind_t v2, orbind_t v3, orbind_t v4, orbind_t v5) const;
  Tensor6_entry_t get(size_t idx) const;
  void set(size_t idx, Tensor6_entry_t val);
  void set(size_t idx, double val);
  double norm() const;
  void print(double thr=1e-10) const;
};
typedef struct {
  orbind_t o0, o1, o2, o3, o4, o5, o6;
  orbind_t v0, v1, v2, v3, v4, v5, v6;
  double el;
} Tensor7_entry_t;
bool operator<(const Tensor7_entry_t & lh, const Tensor7_entry_t & rh);
bool operator==(const Tensor7_entry_t & lh, const Tensor7_entry_t & rh);
Tensor7_entry_t make_entry7(orbind_t o0, orbind_t o1, orbind_t o2, orbind_t o3, orbind_t o4, orbind_t o5, orbind_t o6, orbind_t v0, orbind_t v1, orbind_t v2, orbind_t v3, orbind_t v4, orbind_t v5, orbind_t v6);
void unroll_entry7(const Tensor7_entry_t & etr, orbind_t & o0, orbind_t & o1, orbind_t & o2, orbind_t & o3, orbind_t & o4, orbind_t & o5, orbind_t & o6, orbind_t & v0, orbind_t & v1, orbind_t & v2, orbind_t & v3, orbind_t & v4, orbind_t & v5, orbind_t & v6, double & val);
void unroll_entry7(const Tensor7_entry_t & etr, std::vector<orbind_t> & idx, double & val);
class Tensor7 {
  std::vector<Tensor7_entry_t> list;
public:
  Tensor7();
  Tensor7(const std::vector<determinant_t> & list);
  ~Tensor7();
  void add(const Tensor7_entry_t & key);
  size_t size() const;
  double operator()(orbind_t o0, orbind_t o1, orbind_t o2, orbind_t o3, orbind_t o4, orbind_t o5, orbind_t o6, orbind_t v0, orbind_t v1, orbind_t v2, orbind_t v3, orbind_t v4, orbind_t v5, orbind_t v6) const;
  Tensor7_entry_t get(size_t idx) const;
  void set(size_t idx, Tensor7_entry_t val);
  void set(size_t idx, double val);
  double norm() const;
  void print(double thr=1e-10) const;
};
typedef struct {
  orbind_t o0, o1, o2, o3, o4, o5, o6, o7;
  orbind_t v0, v1, v2, v3, v4, v5, v6, v7;
  double el;
} Tensor8_entry_t;
bool operator<(const Tensor8_entry_t & lh, const Tensor8_entry_t & rh);
bool operator==(const Tensor8_entry_t & lh, const Tensor8_entry_t & rh);
Tensor8_entry_t make_entry8(orbind_t o0, orbind_t o1, orbind_t o2, orbind_t o3, orbind_t o4, orbind_t o5, orbind_t o6, orbind_t o7, orbind_t v0, orbind_t v1, orbind_t v2, orbind_t v3, orbind_t v4, orbind_t v5, orbind_t v6, orbind_t v7);
void unroll_entry8(const Tensor8_entry_t & etr, orbind_t & o0, orbind_t & o1, orbind_t & o2, orbind_t & o3, orbind_t & o4, orbind_t & o5, orbind_t & o6, orbind_t & o7, orbind_t & v0, orbind_t & v1, orbind_t & v2, orbind_t & v3, orbind_t & v4, orbind_t & v5, orbind_t & v6, orbind_t & v7, double & val);
void unroll_entry8(const Tensor8_entry_t & etr, std::vector<orbind_t> & idx, double & val);
class Tensor8 {
  std::vector<Tensor8_entry_t> list;
public:
  Tensor8();
  Tensor8(const std::vector<determinant_t> & list);
  ~Tensor8();
  void add(const Tensor8_entry_t & key);
  size_t size() const;
  double operator()(orbind_t o0, orbind_t o1, orbind_t o2, orbind_t o3, orbind_t o4, orbind_t o5, orbind_t o6, orbind_t o7, orbind_t v0, orbind_t v1, orbind_t v2, orbind_t v3, orbind_t v4, orbind_t v5, orbind_t v6, orbind_t v7) const;
  Tensor8_entry_t get(size_t idx) const;
  void set(size_t idx, Tensor8_entry_t val);
  void set(size_t idx, double val);
  double norm() const;
  void print(double thr=1e-10) const;
};
typedef struct {
  orbind_t o0, o1, o2, o3, o4, o5, o6, o7, o8;
  orbind_t v0, v1, v2, v3, v4, v5, v6, v7, v8;
  double el;
} Tensor9_entry_t;
bool operator<(const Tensor9_entry_t & lh, const Tensor9_entry_t & rh);
bool operator==(const Tensor9_entry_t & lh, const Tensor9_entry_t & rh);
Tensor9_entry_t make_entry9(orbind_t o0, orbind_t o1, orbind_t o2, orbind_t o3, orbind_t o4, orbind_t o5, orbind_t o6, orbind_t o7, orbind_t o8, orbind_t v0, orbind_t v1, orbind_t v2, orbind_t v3, orbind_t v4, orbind_t v5, orbind_t v6, orbind_t v7, orbind_t v8);
void unroll_entry9(const Tensor9_entry_t & etr, orbind_t & o0, orbind_t & o1, orbind_t & o2, orbind_t & o3, orbind_t & o4, orbind_t & o5, orbind_t & o6, orbind_t & o7, orbind_t & o8, orbind_t & v0, orbind_t & v1, orbind_t & v2, orbind_t & v3, orbind_t & v4, orbind_t & v5, orbind_t & v6, orbind_t & v7, orbind_t & v8, double & val);
void unroll_entry9(const Tensor9_entry_t & etr, std::vector<orbind_t> & idx, double & val);
class Tensor9 {
  std::vector<Tensor9_entry_t> list;
public:
  Tensor9();
  Tensor9(const std::vector<determinant_t> & list);
  ~Tensor9();
  void add(const Tensor9_entry_t & key);
  size_t size() const;
  double operator()(orbind_t o0, orbind_t o1, orbind_t o2, orbind_t o3, orbind_t o4, orbind_t o5, orbind_t o6, orbind_t o7, orbind_t o8, orbind_t v0, orbind_t v1, orbind_t v2, orbind_t v3, orbind_t v4, orbind_t v5, orbind_t v6, orbind_t v7, orbind_t v8) const;
  Tensor9_entry_t get(size_t idx) const;
  void set(size_t idx, Tensor9_entry_t val);
  void set(size_t idx, double val);
  double norm() const;
  void print(double thr=1e-10) const;
};
typedef struct {
  orbind_t o0, o1, o2, o3, o4, o5, o6, o7, o8, o9;
  orbind_t v0, v1, v2, v3, v4, v5, v6, v7, v8, v9;
  double el;
} Tensor10_entry_t;
bool operator<(const Tensor10_entry_t & lh, const Tensor10_entry_t & rh);
bool operator==(const Tensor10_entry_t & lh, const Tensor10_entry_t & rh);
Tensor10_entry_t make_entry10(orbind_t o0, orbind_t o1, orbind_t o2, orbind_t o3, orbind_t o4, orbind_t o5, orbind_t o6, orbind_t o7, orbind_t o8, orbind_t o9, orbind_t v0, orbind_t v1, orbind_t v2, orbind_t v3, orbind_t v4, orbind_t v5, orbind_t v6, orbind_t v7, orbind_t v8, orbind_t v9);
void unroll_entry10(const Tensor10_entry_t & etr, orbind_t & o0, orbind_t & o1, orbind_t & o2, orbind_t & o3, orbind_t & o4, orbind_t & o5, orbind_t & o6, orbind_t & o7, orbind_t & o8, orbind_t & o9, orbind_t & v0, orbind_t & v1, orbind_t & v2, orbind_t & v3, orbind_t & v4, orbind_t & v5, orbind_t & v6, orbind_t & v7, orbind_t & v8, orbind_t & v9, double & val);
void unroll_entry10(const Tensor10_entry_t & etr, std::vector<orbind_t> & idx, double & val);
class Tensor10 {
  std::vector<Tensor10_entry_t> list;
public:
  Tensor10();
  Tensor10(const std::vector<determinant_t> & list);
  ~Tensor10();
  void add(const Tensor10_entry_t & key);
  size_t size() const;
  double operator()(orbind_t o0, orbind_t o1, orbind_t o2, orbind_t o3, orbind_t o4, orbind_t o5, orbind_t o6, orbind_t o7, orbind_t o8, orbind_t o9, orbind_t v0, orbind_t v1, orbind_t v2, orbind_t v3, orbind_t v4, orbind_t v5, orbind_t v6, orbind_t v7, orbind_t v8, orbind_t v9) const;
  Tensor10_entry_t get(size_t idx) const;
  void set(size_t idx, Tensor10_entry_t val);
  void set(size_t idx, double val);
  double norm() const;
  void print(double thr=1e-10) const;
};
typedef struct {
  orbind_t o0, o1, o2, o3, o4, o5, o6, o7, o8, o9, o10;
  orbind_t v0, v1, v2, v3, v4, v5, v6, v7, v8, v9, v10;
  double el;
} Tensor11_entry_t;
bool operator<(const Tensor11_entry_t & lh, const Tensor11_entry_t & rh);
bool operator==(const Tensor11_entry_t & lh, const Tensor11_entry_t & rh);
Tensor11_entry_t make_entry11(orbind_t o0, orbind_t o1, orbind_t o2, orbind_t o3, orbind_t o4, orbind_t o5, orbind_t o6, orbind_t o7, orbind_t o8, orbind_t o9, orbind_t o10, orbind_t v0, orbind_t v1, orbind_t v2, orbind_t v3, orbind_t v4, orbind_t v5, orbind_t v6, orbind_t v7, orbind_t v8, orbind_t v9, orbind_t v10);
void unroll_entry11(const Tensor11_entry_t & etr, orbind_t & o0, orbind_t & o1, orbind_t & o2, orbind_t & o3, orbind_t & o4, orbind_t & o5, orbind_t & o6, orbind_t & o7, orbind_t & o8, orbind_t & o9, orbind_t & o10, orbind_t & v0, orbind_t & v1, orbind_t & v2, orbind_t & v3, orbind_t & v4, orbind_t & v5, orbind_t & v6, orbind_t & v7, orbind_t & v8, orbind_t & v9, orbind_t & v10, double & val);
void unroll_entry11(const Tensor11_entry_t & etr, std::vector<orbind_t> & idx, double & val);
class Tensor11 {
  std::vector<Tensor11_entry_t> list;
public:
  Tensor11();
  Tensor11(const std::vector<determinant_t> & list);
  ~Tensor11();
  void add(const Tensor11_entry_t & key);
  size_t size() const;
  double operator()(orbind_t o0, orbind_t o1, orbind_t o2, orbind_t o3, orbind_t o4, orbind_t o5, orbind_t o6, orbind_t o7, orbind_t o8, orbind_t o9, orbind_t o10, orbind_t v0, orbind_t v1, orbind_t v2, orbind_t v3, orbind_t v4, orbind_t v5, orbind_t v6, orbind_t v7, orbind_t v8, orbind_t v9, orbind_t v10) const;
  Tensor11_entry_t get(size_t idx) const;
  void set(size_t idx, Tensor11_entry_t val);
  void set(size_t idx, double val);
  double norm() const;
  void print(double thr=1e-10) const;
};
typedef struct {
  orbind_t o0, o1, o2, o3, o4, o5, o6, o7, o8, o9, o10, o11;
  orbind_t v0, v1, v2, v3, v4, v5, v6, v7, v8, v9, v10, v11;
  double el;
} Tensor12_entry_t;
bool operator<(const Tensor12_entry_t & lh, const Tensor12_entry_t & rh);
bool operator==(const Tensor12_entry_t & lh, const Tensor12_entry_t & rh);
Tensor12_entry_t make_entry12(orbind_t o0, orbind_t o1, orbind_t o2, orbind_t o3, orbind_t o4, orbind_t o5, orbind_t o6, orbind_t o7, orbind_t o8, orbind_t o9, orbind_t o10, orbind_t o11, orbind_t v0, orbind_t v1, orbind_t v2, orbind_t v3, orbind_t v4, orbind_t v5, orbind_t v6, orbind_t v7, orbind_t v8, orbind_t v9, orbind_t v10, orbind_t v11);
void unroll_entry12(const Tensor12_entry_t & etr, orbind_t & o0, orbind_t & o1, orbind_t & o2, orbind_t & o3, orbind_t & o4, orbind_t & o5, orbind_t & o6, orbind_t & o7, orbind_t & o8, orbind_t & o9, orbind_t & o10, orbind_t & o11, orbind_t & v0, orbind_t & v1, orbind_t & v2, orbind_t & v3, orbind_t & v4, orbind_t & v5, orbind_t & v6, orbind_t & v7, orbind_t & v8, orbind_t & v9, orbind_t & v10, orbind_t & v11, double & val);
void unroll_entry12(const Tensor12_entry_t & etr, std::vector<orbind_t> & idx, double & val);
class Tensor12 {
  std::vector<Tensor12_entry_t> list;
public:
  Tensor12();
  Tensor12(const std::vector<determinant_t> & list);
  ~Tensor12();
  void add(const Tensor12_entry_t & key);
  size_t size() const;
  double operator()(orbind_t o0, orbind_t o1, orbind_t o2, orbind_t o3, orbind_t o4, orbind_t o5, orbind_t o6, orbind_t o7, orbind_t o8, orbind_t o9, orbind_t o10, orbind_t o11, orbind_t v0, orbind_t v1, orbind_t v2, orbind_t v3, orbind_t v4, orbind_t v5, orbind_t v6, orbind_t v7, orbind_t v8, orbind_t v9, orbind_t v10, orbind_t v11) const;
  Tensor12_entry_t get(size_t idx) const;
  void set(size_t idx, Tensor12_entry_t val);
  void set(size_t idx, double val);
  double norm() const;
  void print(double thr=1e-10) const;
};
typedef struct {
  orbind_t o0, o1, o2, o3, o4, o5, o6, o7, o8, o9, o10, o11, o12;
  orbind_t v0, v1, v2, v3, v4, v5, v6, v7, v8, v9, v10, v11, v12;
  double el;
} Tensor13_entry_t;
bool operator<(const Tensor13_entry_t & lh, const Tensor13_entry_t & rh);
bool operator==(const Tensor13_entry_t & lh, const Tensor13_entry_t & rh);
Tensor13_entry_t make_entry13(orbind_t o0, orbind_t o1, orbind_t o2, orbind_t o3, orbind_t o4, orbind_t o5, orbind_t o6, orbind_t o7, orbind_t o8, orbind_t o9, orbind_t o10, orbind_t o11, orbind_t o12, orbind_t v0, orbind_t v1, orbind_t v2, orbind_t v3, orbind_t v4, orbind_t v5, orbind_t v6, orbind_t v7, orbind_t v8, orbind_t v9, orbind_t v10, orbind_t v11, orbind_t v12);
void unroll_entry13(const Tensor13_entry_t & etr, orbind_t & o0, orbind_t & o1, orbind_t & o2, orbind_t & o3, orbind_t & o4, orbind_t & o5, orbind_t & o6, orbind_t & o7, orbind_t & o8, orbind_t & o9, orbind_t & o10, orbind_t & o11, orbind_t & o12, orbind_t & v0, orbind_t & v1, orbind_t & v2, orbind_t & v3, orbind_t & v4, orbind_t & v5, orbind_t & v6, orbind_t & v7, orbind_t & v8, orbind_t & v9, orbind_t & v10, orbind_t & v11, orbind_t & v12, double & val);
void unroll_entry13(const Tensor13_entry_t & etr, std::vector<orbind_t> & idx, double & val);
class Tensor13 {
  std::vector<Tensor13_entry_t> list;
public:
  Tensor13();
  Tensor13(const std::vector<determinant_t> & list);
  ~Tensor13();
  void add(const Tensor13_entry_t & key);
  size_t size() const;
  double operator()(orbind_t o0, orbind_t o1, orbind_t o2, orbind_t o3, orbind_t o4, orbind_t o5, orbind_t o6, orbind_t o7, orbind_t o8, orbind_t o9, orbind_t o10, orbind_t o11, orbind_t o12, orbind_t v0, orbind_t v1, orbind_t v2, orbind_t v3, orbind_t v4, orbind_t v5, orbind_t v6, orbind_t v7, orbind_t v8, orbind_t v9, orbind_t v10, orbind_t v11, orbind_t v12) const;
  Tensor13_entry_t get(size_t idx) const;
  void set(size_t idx, Tensor13_entry_t val);
  void set(size_t idx, double val);
  double norm() const;
  void print(double thr=1e-10) const;
};
typedef struct {
  orbind_t o0, o1, o2, o3, o4, o5, o6, o7, o8, o9, o10, o11, o12, o13;
  orbind_t v0, v1, v2, v3, v4, v5, v6, v7, v8, v9, v10, v11, v12, v13;
  double el;
} Tensor14_entry_t;
bool operator<(const Tensor14_entry_t & lh, const Tensor14_entry_t & rh);
bool operator==(const Tensor14_entry_t & lh, const Tensor14_entry_t & rh);
Tensor14_entry_t make_entry14(orbind_t o0, orbind_t o1, orbind_t o2, orbind_t o3, orbind_t o4, orbind_t o5, orbind_t o6, orbind_t o7, orbind_t o8, orbind_t o9, orbind_t o10, orbind_t o11, orbind_t o12, orbind_t o13, orbind_t v0, orbind_t v1, orbind_t v2, orbind_t v3, orbind_t v4, orbind_t v5, orbind_t v6, orbind_t v7, orbind_t v8, orbind_t v9, orbind_t v10, orbind_t v11, orbind_t v12, orbind_t v13);
void unroll_entry14(const Tensor14_entry_t & etr, orbind_t & o0, orbind_t & o1, orbind_t & o2, orbind_t & o3, orbind_t & o4, orbind_t & o5, orbind_t & o6, orbind_t & o7, orbind_t & o8, orbind_t & o9, orbind_t & o10, orbind_t & o11, orbind_t & o12, orbind_t & o13, orbind_t & v0, orbind_t & v1, orbind_t & v2, orbind_t & v3, orbind_t & v4, orbind_t & v5, orbind_t & v6, orbind_t & v7, orbind_t & v8, orbind_t & v9, orbind_t & v10, orbind_t & v11, orbind_t & v12, orbind_t & v13, double & val);
void unroll_entry14(const Tensor14_entry_t & etr, std::vector<orbind_t> & idx, double & val);
class Tensor14 {
  std::vector<Tensor14_entry_t> list;
public:
  Tensor14();
  Tensor14(const std::vector<determinant_t> & list);
  ~Tensor14();
  void add(const Tensor14_entry_t & key);
  size_t size() const;
  double operator()(orbind_t o0, orbind_t o1, orbind_t o2, orbind_t o3, orbind_t o4, orbind_t o5, orbind_t o6, orbind_t o7, orbind_t o8, orbind_t o9, orbind_t o10, orbind_t o11, orbind_t o12, orbind_t o13, orbind_t v0, orbind_t v1, orbind_t v2, orbind_t v3, orbind_t v4, orbind_t v5, orbind_t v6, orbind_t v7, orbind_t v8, orbind_t v9, orbind_t v10, orbind_t v11, orbind_t v12, orbind_t v13) const;
  Tensor14_entry_t get(size_t idx) const;
  void set(size_t idx, Tensor14_entry_t val);
  void set(size_t idx, double val);
  double norm() const;
  void print(double thr=1e-10) const;
};
typedef struct {
  orbind_t o0, o1, o2, o3, o4, o5, o6, o7, o8, o9, o10, o11, o12, o13, o14;
  orbind_t v0, v1, v2, v3, v4, v5, v6, v7, v8, v9, v10, v11, v12, v13, v14;
  double el;
} Tensor15_entry_t;
bool operator<(const Tensor15_entry_t & lh, const Tensor15_entry_t & rh);
bool operator==(const Tensor15_entry_t & lh, const Tensor15_entry_t & rh);
Tensor15_entry_t make_entry15(orbind_t o0, orbind_t o1, orbind_t o2, orbind_t o3, orbind_t o4, orbind_t o5, orbind_t o6, orbind_t o7, orbind_t o8, orbind_t o9, orbind_t o10, orbind_t o11, orbind_t o12, orbind_t o13, orbind_t o14, orbind_t v0, orbind_t v1, orbind_t v2, orbind_t v3, orbind_t v4, orbind_t v5, orbind_t v6, orbind_t v7, orbind_t v8, orbind_t v9, orbind_t v10, orbind_t v11, orbind_t v12, orbind_t v13, orbind_t v14);
void unroll_entry15(const Tensor15_entry_t & etr, orbind_t & o0, orbind_t & o1, orbind_t & o2, orbind_t & o3, orbind_t & o4, orbind_t & o5, orbind_t & o6, orbind_t & o7, orbind_t & o8, orbind_t & o9, orbind_t & o10, orbind_t & o11, orbind_t & o12, orbind_t & o13, orbind_t & o14, orbind_t & v0, orbind_t & v1, orbind_t & v2, orbind_t & v3, orbind_t & v4, orbind_t & v5, orbind_t & v6, orbind_t & v7, orbind_t & v8, orbind_t & v9, orbind_t & v10, orbind_t & v11, orbind_t & v12, orbind_t & v13, orbind_t & v14, double & val);
void unroll_entry15(const Tensor15_entry_t & etr, std::vector<orbind_t> & idx, double & val);
class Tensor15 {
  std::vector<Tensor15_entry_t> list;
public:
  Tensor15();
  Tensor15(const std::vector<determinant_t> & list);
  ~Tensor15();
  void add(const Tensor15_entry_t & key);
  size_t size() const;
  double operator()(orbind_t o0, orbind_t o1, orbind_t o2, orbind_t o3, orbind_t o4, orbind_t o5, orbind_t o6, orbind_t o7, orbind_t o8, orbind_t o9, orbind_t o10, orbind_t o11, orbind_t o12, orbind_t o13, orbind_t o14, orbind_t v0, orbind_t v1, orbind_t v2, orbind_t v3, orbind_t v4, orbind_t v5, orbind_t v6, orbind_t v7, orbind_t v8, orbind_t v9, orbind_t v10, orbind_t v11, orbind_t v12, orbind_t v13, orbind_t v14) const;
  Tensor15_entry_t get(size_t idx) const;
  void set(size_t idx, Tensor15_entry_t val);
  void set(size_t idx, double val);
  double norm() const;
  void print(double thr=1e-10) const;
};
typedef struct {
  orbind_t o0, o1, o2, o3, o4, o5, o6, o7, o8, o9, o10, o11, o12, o13, o14, o15;
  orbind_t v0, v1, v2, v3, v4, v5, v6, v7, v8, v9, v10, v11, v12, v13, v14, v15;
  double el;
} Tensor16_entry_t;
bool operator<(const Tensor16_entry_t & lh, const Tensor16_entry_t & rh);
bool operator==(const Tensor16_entry_t & lh, const Tensor16_entry_t & rh);
Tensor16_entry_t make_entry16(orbind_t o0, orbind_t o1, orbind_t o2, orbind_t o3, orbind_t o4, orbind_t o5, orbind_t o6, orbind_t o7, orbind_t o8, orbind_t o9, orbind_t o10, orbind_t o11, orbind_t o12, orbind_t o13, orbind_t o14, orbind_t o15, orbind_t v0, orbind_t v1, orbind_t v2, orbind_t v3, orbind_t v4, orbind_t v5, orbind_t v6, orbind_t v7, orbind_t v8, orbind_t v9, orbind_t v10, orbind_t v11, orbind_t v12, orbind_t v13, orbind_t v14, orbind_t v15);
void unroll_entry16(const Tensor16_entry_t & etr, orbind_t & o0, orbind_t & o1, orbind_t & o2, orbind_t & o3, orbind_t & o4, orbind_t & o5, orbind_t & o6, orbind_t & o7, orbind_t & o8, orbind_t & o9, orbind_t & o10, orbind_t & o11, orbind_t & o12, orbind_t & o13, orbind_t & o14, orbind_t & o15, orbind_t & v0, orbind_t & v1, orbind_t & v2, orbind_t & v3, orbind_t & v4, orbind_t & v5, orbind_t & v6, orbind_t & v7, orbind_t & v8, orbind_t & v9, orbind_t & v10, orbind_t & v11, orbind_t & v12, orbind_t & v13, orbind_t & v14, orbind_t & v15, double & val);
void unroll_entry16(const Tensor16_entry_t & etr, std::vector<orbind_t> & idx, double & val);
class Tensor16 {
  std::vector<Tensor16_entry_t> list;
public:
  Tensor16();
  Tensor16(const std::vector<determinant_t> & list);
  ~Tensor16();
  void add(const Tensor16_entry_t & key);
  size_t size() const;
  double operator()(orbind_t o0, orbind_t o1, orbind_t o2, orbind_t o3, orbind_t o4, orbind_t o5, orbind_t o6, orbind_t o7, orbind_t o8, orbind_t o9, orbind_t o10, orbind_t o11, orbind_t o12, orbind_t o13, orbind_t o14, orbind_t o15, orbind_t v0, orbind_t v1, orbind_t v2, orbind_t v3, orbind_t v4, orbind_t v5, orbind_t v6, orbind_t v7, orbind_t v8, orbind_t v9, orbind_t v10, orbind_t v11, orbind_t v12, orbind_t v13, orbind_t v14, orbind_t v15) const;
  Tensor16_entry_t get(size_t idx) const;
  void set(size_t idx, Tensor16_entry_t val);
  void set(size_t idx, double val);
  double norm() const;
  void print(double thr=1e-10) const;
};
#endif