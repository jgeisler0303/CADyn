#ifndef EKF_autotune_HPP_
#define EKF_autotune_HPP_

#include "NewmarkBeta.hpp"

class EKFException: public std::exception {
public:
    EKFException(const std::string& msg= "General EKF exception") :
        m_msg(msg)
    {  }
    
    virtual const char* what() const throw () {
        return m_msg.c_str();
    }
    
    const std::string m_msg;
};


template <int nbrstates_, class system_type>
class EKF_autotune {
public:
    EKF_autotune();
    virtual ~EKF_autotune() {}

    static constexpr int nbrdof= system_type::nbrdof;
    static constexpr int nbrin= system_type::nbrin;
    static constexpr int nbrout= system_type::nbrout;
    static constexpr int nbrstates= nbrstates_;

    typedef typename system_type::real_type real_type;
    typedef Eigen::Matrix<real_type, nbrstates_, 1> VecSt;
    typedef typename system_type::VecO VecO;
    typedef Eigen::Matrix<real_type, nbrstates_, nbrstates_> MatSt;
    typedef Eigen::Matrix<real_type, nbrout, nbrout> MatO;
    typedef Eigen::Matrix<real_type, nbrstates_, nbrout> MatStO;
    typedef Eigen::Matrix<real_type, nbrout, nbrstates_> MatOSt;
    typedef Eigen::Matrix<int, nbrdof, 1> VecQ_int;
    
    void next(real_type ts, const Eigen::Ref<const VecO> &y_meas);
    
    system_type system;
    
    
    VecSt x, x_ul, x_ll;
    VecQ_int qx_idx;
    VecQ_int dqx_idx;
    
    MatSt ekfSigma;
    MatSt ekfQ;
    MatO ekfR;
    
    VecO adaptScale;
    VecSt adaptUpdate;
    
    real_type T_adapt= -1.0;
};

// template <int nbrstates_, class system_type> const int EKF_autotune<nbrstates_, system_type>::nbrstates= nbrstates_;
// template <int nbrstates_, class system_type> const int EKF_autotune<nbrstates_, system_type>::nbrdof= system_type::nbrdof;
// template <int nbrstates_, class system_type> const int EKF_autotune<nbrstates_, system_type>::nbrin= system_type::nbrin;
// template <int nbrstates_, class system_type> const int EKF_autotune<nbrstates_, system_type>::nbrout= system_type::nbrout;


template <int nbrstates_, class system_type>
EKF_autotune<nbrstates_, system_type>::EKF_autotune() {
    ekfSigma.setIdentity();
    ekfQ.setIdentity();
    ekfQ*= 1.0e-6;
    ekfR.setIdentity();
    adaptScale.setConstant(1.0);
    adaptUpdate.setConstant(1.0);
    
    x_ul.setConstant(std::numeric_limits<real_type>::infinity());
    x_ll.setConstant(-std::numeric_limits<real_type>::infinity());
}

template <int nbrstates_, class system_type>
void EKF_autotune<nbrstates_, system_type>::next(real_type ts, const Eigen::Ref<const VecO> &y_meas) {
    system.newmarkIntervalWithSens(ts);
    
    MatSt ekfA;
    MatOSt ekfC;

    // compose state vector and Jacobians
    for(int i= 0; i<nbrdof; ++i) {
        if(qx_idx(i)<0) continue;
        
        x(qx_idx(i))= system.q(i);
        
        for(int j= 0; j<nbrdof; ++j)
            if(qx_idx(j)>=0) ekfA(qx_idx(i), qx_idx(j))= system.S(i, j);
        for(int j= 0; j<nbrdof; ++j)
            if(dqx_idx(j)>=0) ekfA(qx_idx(i), dqx_idx(j))= system.S(i, j+nbrdof);
            
        for(int j= 0; j<nbrout; ++j)
            ekfC(j, qx_idx(i))= system.CD(j, i);
    }        
    
    for(int i= 0; i<nbrdof; ++i) {
        if(dqx_idx(i)<0) continue;
        
        x(dqx_idx(i))= system.qd(i);

        for(int j= 0; j<nbrdof; ++j)
            if(qx_idx(j)>=0) ekfA(dqx_idx(i), qx_idx(j))= system.S(i+nbrdof, j);
        for(int j= 0; j<nbrdof; ++j)
            if(dqx_idx(j)>=0) ekfA(dqx_idx(i), dqx_idx(j))= system.S(i+nbrdof, j+nbrdof);

        for(int j= 0; j<nbrout; ++j)
            ekfC(j, dqx_idx(i))= system.CD(j, i+nbrdof);        
    }
    
    // Kalman equation according to Dan Simon eq (7.14)
    MatSt Sigma_pred= ekfA * ekfSigma * ekfA.transpose() + ekfQ;
    MatStO SC= Sigma_pred*ekfC.transpose();
    MatO CSC= ekfC*SC;
    MatO R_= CSC + ekfR;
    MatO R_inv= R_.ldlt().solve(MatO::Identity());
    MatStO ekfK= SC * R_inv;
    
    ekfSigma= Sigma_pred - ekfK*SC.transpose();
    ekfSigma= 0.5*(ekfSigma+ekfSigma.transpose());
    
    VecO d= y_meas-system.y;
    VecSt Dx= ekfK*d;
    
    if(!is_finite(Dx))
        throw EKFException("EKF: Kalman update is not finite.");
        
    x+= Dx.cwiseProduct(adaptUpdate);
    
    // apply limits
    for(int i= 0; i<nbrstates; ++i)
        x(i)= std::min(std::max(x(i), x_ll(i)), x_ul(i));

    // write state vector back to q and dq
    for(int i= 0; i<nbrdof; ++i)
        if(qx_idx(i)>=0) system.q(i)= x(qx_idx(i));
    
    for(int i= 0; i<nbrdof; ++i)
        if(dqx_idx(i)>=0) system.qd(i)= x(dqx_idx(i));
        
    // perform tuning of Q and R
    if(T_adapt>0.0) {
        real_type alpha_adapt= exp(-ts/T_adapt);
        VecO ep= d - ekfC*ekfK*d; // should be real residual: y_meas-h(x_corr), but d= y_meas-y_pred, x_corr=x_pred+K*d, C*K*d=C*x_corr-C*x_pred=y_corr-y_pred, d-C*K*d=y_meas-y_pred-(y_corr-y_pred)=y_meas-y_corr
        ep= ep.cwiseProduct(adaptScale);

        ekfR= alpha_adapt*ekfR + (1.0-alpha_adapt)*(ep*ep.transpose() + CSC);

        ekfQ= alpha_adapt*ekfQ + (1.0-alpha_adapt) * Dx*Dx.transpose();
    }    
}

#endif /* EKF_autotune_HPP_ */
