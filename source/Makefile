# Written by Susi Lehtola, 2016
# Copyright (c) 2016, Susi Lehtola
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.

CXX=g++
CXXFLAGS=-O2 -std=c++11 -g -Wall -fopenmp $(INCLUDE)
GMPLIBS=-lgmpxx -lgmp

all: main_bit.x main_ubit.x generate.x convert_bit_ubit.x convert_ubit_bit.x
all_gmp: $(all) main_int.x convert_bit_int.x convert_ubit_int.x convert_int_bit.x convert_int_ubit.x

tensor.hh tensor.cc: gentensor.py
	python gentensor.py

SOURCES=parser.cc sort.cc tensor.cc wick.cc
OBJS=$(SOURCES:.cc=.o)

#main.cc: generate.x
#	./generate.x

main.o: main.cc evaluate.cc.in gather.hh
main_bit.x: $(OBJS) main.o read_bitstring.o
	$(CXX) $(CXXFLAGS) $(OBJS) main.o read_bitstring.o -o main_bit.x
main_ubit.x: $(OBJS) main.o read_ubitstring.o
	$(CXX) $(CXXFLAGS) $(OBJS) main.o read_ubitstring.o -o main_ubit.x
main_int.x: $(OBJS) main.o read_integerlist.o
	$(CXX) $(CXXFLAGS) $(OBJS) main.o read_integerlist.o $(GMPLIBS) -o main_int.x

generate.x: generate.o wick.o
	$(CXX) $(CXXFLAGS) generate.o wick.o -o generate.x

convert_bit_int.x: convert.o read_bitstring.o write_integerlist.o
	$(CXX) $(CXXFLAGS) $(OBJS) convert.o read_bitstring.o write_integerlist.o $(GMPLIBS) -o convert_bit_int.x
convert_ubit_int.x: convert.o read_ubitstring.o write_integerlist.o
	$(CXX) $(CXXFLAGS) $(OBJS) convert.o read_ubitstring.o write_integerlist.o $(GMPLIBS) -o convert_ubit_int.x
convert_int_bit.x: convert.o read_integerlist.o write_bitstring.o
	$(CXX) $(CXXFLAGS) $(OBJS) convert.o read_integerlist.o write_bitstring.o $(GMPLIBS) -o convert_int_bit.x
convert_ubit_bit.x: convert.o read_ubitstring.o write_bitstring.o
	$(CXX) $(CXXFLAGS) $(OBJS) convert.o read_ubitstring.o write_bitstring.o -o convert_ubit_bit.x
convert_int_ubit.x: convert.o read_ubitstring.o write_ubitstring.o
	$(CXX) $(CXXFLAGS) $(OBJS) convert.o read_integerlist.o write_ubitstring.o $(GMPLIBS) -o convert_int_ubit.x
convert_bit_ubit.x: convert.o read_bitstring.o write_ubitstring.o
	$(CXX) $(CXXFLAGS) $(OBJS) convert.o read_bitstring.o write_ubitstring.o -o convert_bit_ubit.x

%.o: %.cc %.hh
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	\rm *.o
