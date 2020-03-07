/*
 * ODEOrder2.hpp
 *
 *  Created on: 10.03.2017
 *      Author: jgeisler
 */

#ifndef ODEORDER2_HPP_
#define ODEORDER2_HPP_

#include <iosfwd>
#include "MBTypes.hpp"

class AbstractIntegratorVisitor;

class ODEOrder2 {
public:
    ODEOrder2(int nbrdof_, int nbrin_= 0, const string &aname= "anonymous_mbs", const string &adesc= "No description");
    virtual ~ODEOrder2() {}

    void writeStateVariablesHeader(std::ostream &OutFile);
    void writeStateVariables(std::ostream &OutFile);
    VecX computeResidualsInt();
    virtual VecX computeResiduals()= 0;
    virtual void calcJacobian(double alphaM, double alphaC, double alphaK, double tol);
    bool staticEquilibrium();
    int newmarkOneStep(double h, bool hmodified= true);
    bool newmarkInterval(double tfinal, double &h, double hmax);
    bool newmarkIntegration(double tfinal, double hsave, double hmax, AbstractIntegratorVisitor *visitor= nullptr);

    string name;
    string description;

    const int nbrdof, nbrin;
    std::vector<string> state_name;
    std::vector<string> in_name;
    MatX Jacobian;
    Eigen::FullPivLU<MatX> LU;
    VecX q, qd, qdd;
    VecX f;
    std::vector<bool> doflocked;
    VecX u;
    double t;
    double AbsTol= 1E-6;
    double RelTol= 1E-6;
    double StepTol= 1E-8;
    double hminmin= 1E-8;
    int jac_recalc_step= 4;
    int max_steps= 10;
    
    int n_back_steps;
    int n_steps;
    double errq;
};

class AbstractIntegratorVisitor {
public:
    AbstractIntegratorVisitor(): system(nullptr) {}
    virtual ~AbstractIntegratorVisitor() {}

    void setSystem(ODEOrder2 *system_) { system= system_; }
    virtual void start()= 0;
    virtual void step()= 0;
    virtual void finish()= 0;

    ODEOrder2 *system;
};

#endif /* ODEORDER2_HPP_ */
