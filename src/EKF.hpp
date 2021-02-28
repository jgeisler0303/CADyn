/*
 *      Author: jgeisler
 */

#ifndef EKF_HPP_
#define EKF_HPP_

#include "MBTypes.hpp"

class EKF {
public:
    EKF(int n_states, int n_out);
    virtual ~EKF() {}

    void next(const Eigen::Ref<const VecX> &u, const Eigen::Ref<const VecX> &y_meas);
    virtual void callModel(const Eigen::Ref<const VecX> &u) = 0;
    
    int n_states;
    int n_out;
    
    VecX x;
    VecX y;
    
    MatX ekfSigma;
    MatX ekfQ;
    MatX ekfR;
    MatX ekfN;
    MatX ekfA;
    MatX ekfC;
};

EKF::EKF(int n_states, int n_out) :
        n_states(n_states),
        n_out(n_out),
        x(n_states),
        y(n_out),
        ekfSigma(n_states, n_states),
        ekfQ(n_states, n_states),
        ekfR(n_out, n_out),
        ekfN(n_states, n_out),
        ekfA(n_states, n_states),
        ekfC(n_out, n_states)
{
    ekfSigma.setIdentity();
}

void EKF::next(const Eigen::Ref<const VecX> &u, const Eigen::Ref<const VecX> &y_meas) {
    // call Model must update x with a prediction, y with the output and possibly ekfA, ekfC, ekfQ, ekfR and ekfN
    callModel(u);
    
    MatX ASigma= ekfA*ekfSigma.selfadjointView<Eigen::Lower>();
    MatX ASigmaCN= ASigma*ekfC.transpose() + ekfN;
    MatX K= ASigmaCN * (ekfC*ekfSigma.selfadjointView<Eigen::Lower>()*ekfC.transpose() + ekfR).selfadjointView<Eigen::Lower>().ldlt().solve(MatX::Identity(n_out, n_out));
    ekfSigma.triangularView<Eigen::Lower>()= ASigma*ekfA.transpose() + ekfQ - K*ASigmaCN.transpose();
    
    x+= K*(y_meas-y);
}

#endif /* EKF_HPP_ */
