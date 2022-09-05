/*
 *      Author: jgeisler
 */

#ifndef EKF_HPP_
#define EKF_HPP_

template <int nbrdof_, int nbrin_, int nbrout_, int nbrstates_, class real_type, class system_type>
class EKF {
public:
    EKF();
    virtual ~EKF() {}

    static const int nbrdof;
    static const int nbrin;
    static const int nbrout;
    static const int nbrstates;

    typedef Eigen::Matrix<real_type, nbrstates_, 1> VecSt;
    typedef typename system_type::VecO VecO;
    typedef Eigen::Matrix<real_type, nbrstates_, nbrstates_> MatSt;
    typedef Eigen::Matrix<real_type, nbrout_, nbrout_> MatO;
    typedef Eigen::Matrix<real_type, nbrstates_, nbrout_> MatStO;
    typedef Eigen::Matrix<real_type, nbrout_, nbrstates_> MatOSt;
    typedef Eigen::Matrix<int, nbrdof_, 1> VecQ_int;
    
    bool next(real_type ts, const Eigen::Ref<const VecO> &y_meas);
    
    system_type system;
    
    
    VecSt x, x_ul, x_ll;
    VecQ_int qx_idx;
    VecQ_int dqx_idx;
    
    MatSt ekfSigma;
    MatSt ekfQ;
    MatO ekfR;
    MatStO ekfN;
};

template <int nbrdof_, int nbrin_, int nbrout_, int nbrstates_, class real_type, class system_type> const int EKF<nbrdof_, nbrin_, nbrout_, nbrstates_, real_type, system_type>::nbrstates= nbrstates_;
template <int nbrdof_, int nbrin_, int nbrout_, int nbrstates_, class real_type, class system_type> const int EKF<nbrdof_, nbrin_, nbrout_, nbrstates_, real_type, system_type>::nbrdof= nbrdof_;
template <int nbrdof_, int nbrin_, int nbrout_, int nbrstates_, class real_type, class system_type> const int EKF<nbrdof_, nbrin_, nbrout_, nbrstates_, real_type, system_type>::nbrin= nbrin_;
template <int nbrdof_, int nbrin_, int nbrout_, int nbrstates_, class real_type, class system_type> const int EKF<nbrdof_, nbrin_, nbrout_, nbrstates_, real_type, system_type>::nbrout= nbrout_;


template <int nbrdof_, int nbrin_, int nbrout_, int nbrstates_, class real_type, class system_type>
EKF<nbrdof_, nbrin_, nbrout_, nbrstates_, real_type, system_type>::EKF() {
    ekfSigma.setIdentity();
    x_ul.setConstant(std::numeric_limits<real_type>::infinity());
    x_ll.setConstant(-std::numeric_limits<real_type>::infinity());
}

template <int nbrdof_, int nbrin_, int nbrout_, int nbrstates_, class real_type, class system_type>
bool EKF<nbrdof_, nbrin_, nbrout_, nbrstates_, real_type, system_type>::next(real_type ts, const Eigen::Ref<const VecO> &y_meas) {
    bool res= system.newmarkIntervalWithSens(ts);
    
    if(!res)
        return res;
    
    MatSt ekfA;
    MatOSt ekfC;

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
    
    // Dan Simon eq (7.14)
    MatSt Sigma_pred= ekfA * ekfSigma * ekfA.transpose() + ekfQ;
    MatStO SC= Sigma_pred*ekfC.transpose();
    MatStO ekfK= SC + ekfN;
    MatO ekfCN= ekfC*ekfN;
    MatO R_= ekfC*SC + ekfCN + ekfCN.transpose() + ekfR;
    MatO R_inv= R_.ldlt().solve(MatO::Identity());
    ekfK*= R_inv;
    
    ekfSigma= Sigma_pred - ekfK*(SC.transpose() + ekfN.transpose());
    ekfSigma= 0.5*(ekfSigma+ekfSigma.transpose());
    
    x+= ekfK*(y_meas-system.y);
    
    for(int i= 0; i<nbrstates; ++i)
        x(i)= std::min(std::max(x(i), x_ll(i)), x_ul(i));

    for(int i= 0; i<nbrdof; ++i)
        if(qx_idx(i)>=0) system.q(i)= x(qx_idx(i));
    
    for(int i= 0; i<nbrdof; ++i)
        if(dqx_idx(i)>=0) system.qd(i)= x(dqx_idx(i));
        
    return true;
}

#endif /* EKF_HPP_ */
