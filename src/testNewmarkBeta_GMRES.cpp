#include "NewmarkBeta_GMRES.hpp"
#include <iostream>

class mySolver : public NewmarkBeta_GMRES<2, 1, 1, double> {
public:    
    virtual VecQ computeResiduals(const VecQ &qdd, const VecQ &qd, const VecQ &q) const;
    virtual VecQ calcJacobianProd(const VecQ &b) const;
    virtual void calcOut();

    
    const double m1= 1.0;
    const double m2= 2.0;
    const double k1= 3.0;
    const double k2= 10.0;
    const double k12= 6.0;
    const double d1= 0.5;
    const double d2= 0.8;
};

mySolver::VecQ mySolver::computeResiduals(const VecQ &qdd, const VecQ &qd, const VecQ &q) const {
    VecQ f;
    double qd0= pow(fabs(qd[0]), 3.0);
    if(qd[0]<0.0) qd0= -qd0;
    
    f[0]= m1*qdd[0] + d1*qd0 + k1*q[0] + k12*(q[0]-q[1]);
    f[1]= m2*qdd[1] + d2*qd[1] + k2*q[1] + k12*(q[1]-q[0]);
    
    return f;
}

mySolver::VecQ mySolver::calcJacobianProd(const VecQ &b) const {
    VecQ Dq= alphaK * b;
    VecQ Dqd= alphaC * b;
    VecQ Dqdd= alphaM * b;

    VecQ df;
    
    df[0]= m1*Dqdd[0] + 3.0*d1*pow(fabs(qd[0]), 2.0)*Dqd[0] + k1*Dq[0] + k12*(Dq[0]-Dq[1]);
    df[1]= m2*Dqdd[1] + d2*Dqd[1] + k2*Dq[1] + k12*(Dq[1]-Dq[0]);
    
    return df;
}

void mySolver::calcOut() {
    y[0]= q[0];
}


mySolver sol;

int main() {
    double h= 0.1;
    sol.q[0]= 1.0;
    
//     std::cout << "q0: " << sol.q << std::endl;
//     std::cout << "qd0: " << sol.qd << std::endl;
//     std::cout << "qdd0: " << sol.qdd << std::endl;
    
    sol.newmarkOneStep(0.0);
//     std::cout << "q0: " << sol.q << std::endl;
//     std::cout << "qd0: " << sol.qd << std::endl;
//     std::cout << "qdd0: " << sol.qdd << std::endl;
//     sol.newmarkInterval(0.1, h, 0.1);
    sol.newmarkInterval(1.0, h, 0.1);
    sol.calcOut();
    
    std::cout << "n_steps: " << sol.n_steps << std::endl;
    std::cout << "t: " << sol.t << std::endl;
    std::cout << "y: " << sol.y[0] << std::endl;
}
