#!/usr/bin/perl

$CXX = 'g++-4.2';
#$SYS = '-isysroot /Developer/SDKs/MacOSX10.5.sdk';
#$MACHOPT = '-DMACOSX_DEPLOYMENT_TARGET=10.5 -mmacosx-version-min=10.5';

$INCLUDE_PATH = '/usr/local/include/openbabel-2.0';
#$LIBRARY_PATH = '/Users/francisgaudreault/Development/NRGsuite/Executables/MacOSX5+/FlexAID/WRK/libs';
$LIBRARY_PATH = '/Users/francisgaudreault/Programs/openbabel-2.3.2/10.5/lib';

$cmd = "$CXX -Wall -o Process_Ligand Process_Ligand.cpp geometry.cpp -g -Wall -I$INCLUDE_PATH -L$LIBRARY_PATH -lopenbabel $SYS -arch x86_64 $MACHOPT";

print "$cmd\n";
`$cmd`;
