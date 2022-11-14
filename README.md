# mcaa

## REQUIREMENTS

- New-ish c++ compiler that supports c++20 features (used mainly for multithreading support)
- cmake
- tbb (libtbb-dev in apt)
- boost library (libboost-dev in apt) 
- GSL (libgsl-dev in apt)
- LHAPDF

## BUILD

Building the simulator should be as simple as running `cmake . && cmake --build .` in the root directory. If your `lhapdf-config` is not working properly (I've had a lot of issues with it at times), you may have to  modify two files manually: `CMakelists.txt` and `src/CMakelists.txt`. Check the comments in those files.

If you are having trouble with the dependencies, consult `fresh_install.txt`, which is a listing of all the commands I used to get the simulator up and running on a freshly installed ubuntu system