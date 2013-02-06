#!/usr/bin/perl

$CXX = 'g++-4.0';
$SYS = '-isysroot /Developer/SDKs/MacOSX10.4u.sdk';
$MACHOPT = '-DMACOSX_DEPLOYMENT_TARGET=10.4 -mmacosx-version-min=10.4';


$INCLUDE_PATH = '/usr/local/include/openbabel-2.0';
$LIBRARY_PATH = '/Users/francisgaudreault/Development/NRGsuite/Executables/MacOSX4+';

`$CXX -Wall -o Process_Ligand Process_Ligand.cpp geometry.cpp -O3 -I$INCLUDE_PATH -L$LIBRARY_PATH -lopenbabel $SYS -arch i386 $MACHOPT`;
