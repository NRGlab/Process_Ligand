#!/usr/bin/perl

$CXX = 'g++-apple-4.2';

$INCLUDE_PATH = '/opt/local/include/openbabel-2.0';
$LIBRARY_PATH = '/opt/local/lib';

`$CXX -Wall -o Process_Ligand Process_Ligand.cpp geometry.cpp -O3 -I$INCLUDE_PATH -L$LIBRARY_PATH -lopenbabel`;
