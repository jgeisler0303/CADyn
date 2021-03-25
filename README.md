Introduction
============
CADyn stands for Computer Aided Dynamics and is a framwork for generating custom-code to simulate multi-body systems dynamic behavior. It is a complete rework of the [EasyDyn](http://hosting.umons.ac.be/html/mecara/EasyDyn/) software.
The purpose of the rework was firstly to provide an object oriented architecture and interface using modern C++ that encapsulates all internal aspects of the simulation process, and secondly to replace the dependency to the GSL library for vector and matrix calculations by the more C++ friendly [eigen3](http://eigen.tuxfamily.org) library.
Although the algorithms used remain basically the same as in EasyDyn, all code was rewritten from scratch. This and the new interface together with some new features were reasons to give this project a new name.

In addition to the reworked C++ code this framework also features a completely new version of the script for generating the custom code. The script that was formerly written for MuPAD was ported to [Maxima](http://maxima.sourceforge.net/) in order to make it usable again with a completely free software tool-chain.

For the simulation of multi-body systems the script cagem (Computer-Aided Generation of Motion) is used to generate problem specific code of system's dynamic equations of motion form an analytic description of the system configuration (position and pose of each body) of the system. The description of the system is written in terms of the minimal coordinates (degrees of freedom) of the system thus eliminating the need to explicitly state (and handle) the constraints of the relative motion between the bodies. The equations of motion are analytically derived via the [principal of virtual work](https://en.wikipedia.org/wiki/Virtual_work) using the symbolic computation and simplification capabilities of the Computer Algebra System (CAS) running the script (in this case now Maxima). This combination of minimal coordinates and symbolically optimized equations generated into custom C++ code yields a very fast and computationally efficient means to simulate the behavior of multi-body systems.

Beyond the former capability of EasyDyn, CADyn also supports a first basic approach to model systems comprising flexible bodies. The structured object oriented architecture allows to extends this capability with more detailed representations of flexible bodies in the future.

Usage
===============
The definition of the system to be simulated is written in the [Maxima language](http://maxima.sourceforge.net/docs/manual/maxima_toc.html#SEC_Contents) using helper functions for the formulation of the homogeneous transformations defining the bodies' positions and poses in dependence of the minimal coordinates  `q` (aka configuration variables or degrees of freedom). These helper functions are defined in the script file `cagem.mac`. Please take a look at the example files to get a feeling how this is done. A more complete manual may be written later. Meanwhile you may also refer to the manual supplied in the [download of EasyDyn](http://hosting.umons.ac.be/html/mecara/EasyDyn/EasyDyn-1.2.4.tgz) or [here](http://hosting.umons.ac.be/html/mecara/EasyDyn/EasyDyn.pdf).

The syntax of MuPAD and Maxima is very similar. The major difference is that definitions in Maxima are indicated by a colon (:) and completed by a semicolon (;). Further, the array index base was changed from zero to one. So, mass, inertia terms (Ixx, Iyy, Izz, Ixy, Ixz, Iyz), body pose (`T0G`) and minimal coordinates (`q`) all start with `mass[1]`, `Ixx[1]`, etc., `T0G[1]` for the first body resp. `q[1]` for the first configuration variable. Nonetheless, the C++ representations of the corresponding entities have a shifted, zero-based index.

The problem-specific code is generated with the function `cagem` also defined in the script file `cagem.mac`. The actual code is written using the the maxima package `gentran`. Unfortunately, the gentran package as supplied by the current maxima release is not fully functional. Therefore you will have to download the package from [my repo](https://github.com/jgeisler0303/maxima) and replace the files in the `maxima/share/contrib/gentran` folder of your maxima installation with the files from my repo.

Further, in the file `maxima/share/numeric/expense.lisp` the line `(declare (fixnum exin rem in-cut pw2))` must be commented out (put a `;` at the start of the line).

The Maxima installation must be compiled for [Steel Bank Common Lisp (SBCL)](http://www.sbcl.org/). On Windows, SBCL is the default Lisp implementation supplied with the [wxMaxima](http://andrejv.github.io/wxmaxima/) installation. For Linux you can download appropriate installation packages from the Sourcforge home page of the Maxima project.

The function `cagem` is called with the file name of the Maxima script defining the system as the first argument and an optional second argument giving the desired path for the generated code. The function generates the main code for the system object in a file with the name of the script defining the system postfixed "System.hpp" and one file postfixed "_app.cpp" that contain an executable test program for the generated system. If present, a file in the same directory as the script defining the system and postfixed ".AppEff.hpp" is merged into the test program. This allows for the definition of time varying external forces on the system or forces between the bodies of the system such as spring or damper forces.

The concrete workflow is as follows:

1. Write the system definition in the maxima script language defining all relevant information expected by the generator script.
1. Start Maxima or wxMaxima.
1. Load the cagem script: `load("cagem.mac");` (assuming you are in the folder where the script resides, otherwise add the path to the argument of the load command).
1. Run the generation script: `cagem("path_to_system_def/{filename_of_system}.mac");`. The custom code will be generated in the current folder with the name `{filename_of_system}System.hpp` and an executable test program is generated with the name `{filename_of_system}_app.cpp`.
1. Compile the custom code and link it to the CADyn library e.g. using gcc: `g++ -std=c++0x -I. -Ipath_to_CADyn/src {filename_of_system}_app.cpp -Lpath_to_CADyn_Library -lCADyn -o {filename_of_system}_app`. This assumes the [eigen3](http://eigen.tuxfamily.org) library is installed and on your system include path.
1. Now you can run the simulation of your system according to the initial conditions and duration specified in the system definition script. This will create an ASCII file `{filename_of_system}.res` that holds the time series of the system configuration and their time derivatives. Additionally the file `{filename_of_system}.plt` is generated which a [GNU Plot](http://www.gnuplot.info/) script to create some neat PostScript figures of the simulated evolution of the system configuration (run `gnuplot {filename_of_system}.plt` to do this).

Alternatively, you can copy the `makefile` from the `examples` directory into the directory where your system description  -- and only this -- lives. There you can run `make all` which will perform all above steps in one.

This instruction assumes that you built the CADyn library using the provided eclipse project or by running `make all` in the "Debug" directory.
