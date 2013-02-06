#!/usr/bin/perl

$CXX = 'g++-4.2';
$SYS = '-isysroot /Developer/SDKs/MacOSX10.5.sdk';
$MACHOPT = '-DMACOSX_DEPLOYMENT_TARGET=10.5 -mmacosx-version-min=10.5';


$INCLUDE_PATH = '/usr/local/include/openbabel-2.0';
$LIBRARY_PATH = '/Users/francisgaudreault/Development/NRGsuite/Executables/MacOSX5+';

$dir = '/Users/francisgaudreault/Development/Atom_Types';

`$CXX -Wall -o Process_Ligand Process_Ligand.cpp geometry.cpp -O3 -I$INCLUDE_PATH -L$LIBRARY_PATH -lopenbabel $SYS -arch i386 $MACHOPT -arch x86_64 $MACHOPT`;

`ln -fs $dir/src/Process_Ligand $dir/Process_Ligand`;
