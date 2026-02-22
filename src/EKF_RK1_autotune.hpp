#ifndef EKF_RK1_autotune_HPP_
#define EKF_RK1_autotune_HPP_

#include "RK1condensed.hpp"

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


template <class problem_class>
class EKF_RK1_autotune {
public:
    static constexpr auto dimI  = problem_class::dimI;
    static constexpr auto dimII = problem_class::dimII;
    static constexpr auto dimX  = problem_class::dimX;
    static constexpr auto dimIn  = problem_class::dimIn;
    static constexpr auto dimOut = problem_class::dimOut;

    using real_type = typename problem_class::real_type;
    using VecI = typename problem_class::VecI;
    using VecII = typename problem_class::VecII;
    using VecX = typename problem_class::VecX;
    using MatI = typename problem_class::MatI;
    using MatII = typename problem_class::MatII;
    using MatII_I = typename problem_class::MatII_I;
    using MatX = typename problem_class::MatX;
    using MatIn = typename problem_class::MatIn;
    using VecIn = typename problem_class::VecIn;
    using MatII_In = typename problem_class::MatII_In;
    using VecOut = typename problem_class::VecOut;
    using MatOutX = typename problem_class::MatOutX;
    using MatOutII = typename problem_class::MatOutII;
    using MatOutIn = typename problem_class::MatOutIn;
    typedef Eigen::Matrix<real_type, dimOut, dimOut> MatOut;
    typedef Eigen::Matrix<real_type, dimX, dimOut> MatXOut;

    EKF_RK1_autotune();
    virtual ~EKF_RK1_autotune() {}

    void next(real_type ts, const Eigen::Ref<const VecOut> &y_meas);
    
    // TODO: consider making this a base class
    RK1condensed<problem_class> rk1_solver;
    
    VecX x_ul, x_ll; // clamping limits for state vector
    
    MatX ekfSigma;
    MatX ekfSigma_pred;
    MatX ekfQ;
    MatOut ekfR;
    MatOut ekfS;
    
#ifdef ZERO_MEAN_TUNING
    VecOut ep_mean;
    VecX Dx_mean;
#endif
    
    real_type d_norm;
    
    VecOut adaptScale;
    VecX fixedQxx;
    VecOut fixedRxx;
    
    real_type T_adapt= -1.0;
};

template <class problem_class>
EKF_RK1_autotune<problem_class>::EKF_RK1_autotune() {
    ekfSigma.setIdentity();
    ekfQ.setIdentity();
    ekfQ*= 1.0e-6;
    ekfR.setIdentity();
    adaptScale.setConstant(1.0);
    fixedQxx.setZero();
    fixedRxx.setZero();
    
    x_ul.setConstant(std::numeric_limits<real_type>::infinity());
    x_ll.setConstant(-std::numeric_limits<real_type>::infinity());

#ifdef ZERO_MEAN_TUNING
    ep_mean.setZero();
    Dx_mean.setZero();
#endif
}

template <class problem_class>
void EKF_RK1_autotune<problem_class>::next(real_type ts, const Eigen::Ref<const VecOut> &y_meas) {
    rk1_solver.intervalWithSens(ts); // this updates rk1_solver.y, rk1_solver.fx, rk1_solver.gx
    
    // Kalman equation according to Dan Simon eq (7.14)
    ekfSigma_pred= rk1_solver.fx * ekfSigma * rk1_solver.fx.transpose() + ekfQ;
    MatXOut SigC= ekfSigma_pred*rk1_solver.gx.transpose();
    MatOut CSigC= rk1_solver.gx*SigC;
    ekfS= CSigC + ekfR;
    Eigen::LDLT<MatOut> ldltS(ekfS);
    if(ldltS.info()!=Eigen::Success)
        throw EKFException("EKF: Innovation covariance decomposition failed.");
    MatXOut ekfK= ldltS.solve(SigC.transpose()).transpose();
    
    // Joseph form for numerical stability: Sigma= (I-KC)*Sigma_pred*(I-KC)' + K*R*K'
    MatX I_KC= MatX::Identity() - ekfK*rk1_solver.gx;
    ekfSigma= I_KC*ekfSigma_pred*I_KC.transpose() + ekfK*ekfR*ekfK.transpose();
    ekfSigma= 0.5*(ekfSigma+ekfSigma.transpose());  // ensure symmetry

    VecOut d= y_meas - rk1_solver.y;
    d_norm= d.dot(ldltS.solve(d));
    
    VecX Dx= ekfK*d;
    
    if(!is_finite(Dx))
        throw EKFException("EKF: Kalman update is not finite.");
        
    rk1_solver.x+= Dx;
    
    // apply limits
    rk1_solver.x = rk1_solver.x.cwiseMax(x_ll).cwiseMin(x_ul);

    // perform tuning of Q and R
    if(T_adapt>0.0) {
        real_type alpha_adapt= exp(-ts/T_adapt);
        VecOut ep= d - rk1_solver.gx*ekfK*d; // should be real residual: y_meas-h(x_corr), but d= y_meas-y_pred, x_corr=x_pred+K*d, C*K*d=C*x_corr-C*x_pred=y_corr-y_pred, d-C*K*d=y_meas-y_pred-(y_corr-y_pred)=y_meas-y_corr
#ifdef ZERO_MEAN_TUNING
        ep_mean= alpha_adapt*ep_mean + (1.0-alpha_adapt)*ep;
        ep-= ep_mean;
#endif        
        ep= ep.cwiseProduct(adaptScale);

        ekfR= alpha_adapt*ekfR + (1.0-alpha_adapt)*(ep*ep.transpose() + CSigC);
        for(int i= 0; i<dimOut; ++i) {
            if(fixedRxx(i)>0.0) {
                ekfR.row(i).setZero();
                ekfR.col(i).setZero();
                ekfR(i, i)= fixedRxx(i);
            }
            else if(fixedRxx(i)<0.0) {
                ekfR(i, i)= -fixedRxx(i);
            }
        }

#ifdef ZERO_MEAN_TUNING
        Dx_mean= alpha_adapt*Dx_mean + (1.0-alpha_adapt)*Dx;        
        Dx-= Dx_mean;
#endif        
        ekfQ= alpha_adapt*ekfQ + (1.0-alpha_adapt) * Dx*Dx.transpose();
        for(int i= 0; i<dimX; ++i) {
            if(fixedQxx(i)>0.0) {
                ekfQ.row(i).setZero();
                ekfQ.col(i).setZero();
                ekfQ(i, i)= fixedQxx(i);
            }
            else if(fixedQxx(i)<0.0) {
                ekfQ(i, i)= -fixedQxx(i);
            }
        }

    }    
}

#endif /* EKF_RK1_autotune_HPP_ */
