Compile ProcessLigand using CMake:
```
git clone --branch cmake https://github.com/NRGlab/Process_Ligand
cd Process_Ligand
mkdir build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release -DBUILD_STATIC_EXECUTABLE=ON
cmake --build . --target ProcessLigand -j 4
```

Replace 4 in the last command with the number of cores to use for compiling in parallel.