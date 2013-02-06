#!/usr/bin/perl

$CXX = 'g++';

$INCLUDE_PATH = '/usr/local/include/openbabel-2.0';
$LIBRARY_PATH = '/usr/local/lib';

`$CXX -Wall -o Process_Ligand Process_Ligand.cpp geometry.cpp -O3 -I$INCLUDE_PATH -L$LIBRARY_PATH -lopenbabel`;
