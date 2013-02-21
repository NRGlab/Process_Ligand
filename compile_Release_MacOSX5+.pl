#!/usr/bin/perl

$CXX = 'g++-4.2';
$SYS = '-isysroot /Developer/SDKs/MacOSX10.5.sdk';
$MACHOPT = '-DMACOSX_DEPLOYMENT_TARGET=10.5 -mmacosx-version-min=10.5';


$INCLUDE_PATH = '/Users/francisgaudreault/Programs/openbabel-2.3.2/10.5/include/openbabel-2.0';
$LIBRARY_PATH = '/Users/francisgaudreault/Programs/openbabel-2.3.2/10.5/lib';

`$CXX -Wall -o Process_Ligand Process_Ligand.cpp geometry.cpp -O3 -I$INCLUDE_PATH -L$LIBRARY_PATH -lopenbabel $SYS -arch i386 $MACHOPT -arch x86_64 $MACHOPT`;
