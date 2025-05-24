# MC-EKRT

MC-EKRT is a Monte-Carlo implementation of the EKRT model for computing partonic initial states in high-energy nuclear collisions. The MC-EKRT event generator is based on collinearly factorized pQCD minijet production, supplemented with a saturation conjecture to control low-pT particle production. It gives a  3-dimensional partonic initial state event-by-event, accounts for dynamical minijet-multiplicity fluctuations in the saturation and particle production, introduces a new type of spatially dependent nuclear parton distribution functions, and accounts for the energy/momentum and valence-quark-number conservations.

When you use the code, please cite the following publications:

[1] M. Kuha, J. Auvinen, K. J. Eskola, H. Hirvonen, Y. Kanakubo and H. Niemi,
"MC-EKRT: Monte Carlo event generator with saturated minijet production for initializing 3+1 D fluid dynamics in high energy nuclear collisions,"
Phys. Rev. C 111, no.5, 054914 (2025)
doi:10.1103/PhysRevC.111.054914
[arXiv:2406.17592 [hep-ph]]

[2] H. Hirvonen, M. Kuha, J. Auvinen, K. J. Eskola, Y. Kanakubo and H. Niemi,
"Effects of saturation and fluctuating hotspots for flow observables in ultrarelativistic heavy-ion collisions,"
Phys. Rev. C 110, no.3, 034911 (2024)
doi:10.1103/PhysRevC.110.034911
[arXiv:2407.01338 [hep-ph]]

More details about the MC-EKRT model and the code can be found in Mikko Kuha's PhD thesis:
M. Kuha, "Simulation of Quark-gluon Plasma Initial States with Monte-Carlo EKRT Model," 2024
https://urn.fi/URN:ISBN:978-952-86-0330-6

The DOI for the version 1.0.0 of the code introduced and used in the above publications is
https://doi.org/10.5281/zenodo.15349744


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
