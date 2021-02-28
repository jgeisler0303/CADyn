#include <iostream>
#include <stdlib.h>
#include <chrono>
#include <cstdio>
#include <ctime>
#include <cmath>
// #include <cxxopts.hpp>

#include "EKF.hpp"

class myEKF: public EKF {
public:
    myEKF() : EKF(2, 1) {
        ekfA(0, 0)= 1.8097;
        ekfA(0, 1)= -0.8187;
        ekfA(1, 0)= 1.0;
        ekfA(1, 1)= 0.0;
        
        ekfC(0, 0)= 0.0374;
        ekfC(0, 1)= 0.0350;
        
        ekfQ.setIdentity();
        ekfR.setIdentity();
        ekfN.setZero();
        
        x.setZero();
    };
    
    virtual void callModel(const Eigen::Ref<const VecX> &u) {
        y= ekfC*x;
        
        x= ekfA*x;
        x(0)+= 0.1250*u(0);
    }
};

int main(int argc, char* argv[]) {
//     cxxopts::Options argc_options("TurbineSimulator", "A simple wind turbine simulator");
//     argc_options.add_options()
//     // ("i,icfile", "Initial conditions file name", cxxopts::value<std::string>()->default_value("./icfile.txt"))
//     ("p,paramfile", "Parameter file name", cxxopts::value<std::string>()->default_value("./params.txt"))
//     ("t,simtime", "Simulation time", cxxopts::value<double>()->default_value("10.0"))
//     ("s,simstep", "Simulation time", cxxopts::value<double>()->default_value("0.01"))
//     ("w,vwind", "Inflow wind definition file name", cxxopts::value<std::string>()->default_value("Inflow.dat"))
//     ("d,discon_dll", "Path and name of the DISCON controller DLL", cxxopts::value<std::string>()->default_value("./discon.dll"))
//     ("o,output", "Output file name", cxxopts::value<std::string>()->default_value("sim_output.outb"))
//     ("a,adjust_wind", "Adjustment factor for wind speed", cxxopts::value<double>()->default_value("1.0"))
//     ;
    
//     auto argc_result = argc_options.parse(argc, argv);
    myEKF ekf;
    
    VecX u(1);
    VecX y(1);
    
    u(0)= 1.0;
    y(0)= 1.0;
    
    for(int i= 0; i<10; ++i) {
        ekf.next(u, y);
        std::cout << "x= " << ekf.x.transpose() << std::endl;
    }
    
    exit (EXIT_SUCCESS);
}
