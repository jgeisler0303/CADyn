/*
 * NewmarkBeta.hpp
 *
 *  Created on: 10.03.2017
 *      Author: jgeisler
 */

#ifndef NEWMARKBETA_HPP_
#define NEWMARKBETA_HPP_

#include <exception>
#include <iosfwd>
#include <cmath>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <vector>
#include <array>
#include <Eigen/Dense>
#include <Eigen/Geometry>

class NewmarkBetaException: public std::exception {
public:
    NewmarkBetaException(const std::string& msg= "General NewmarkBeta integration exception") :
        m_msg(msg)
    {  }
    
    virtual const char* what() const throw () {
        return m_msg.c_str();
    }
    
    const std::string m_msg;
};

template<typename Derived>
inline bool is_finite(const Eigen::MatrixBase<Derived>& x)
{
	return ( (x - x).array() == (x - x).array()).all();
}

template <int nbrdof_, int nbrin_, int nbrout_, class real_type_>
class NewmarkBeta {
public:
    class AbstractIntegratorVisitor {
    public:
        AbstractIntegratorVisitor(): system(nullptr) {}
        virtual ~AbstractIntegratorVisitor() {}

        void setSystem(NewmarkBeta<nbrdof_, nbrin_, nbrout_, real_type_> *system_) { system= system_; }
        virtual void start()= 0;
        virtual void step()= 0;
        virtual void finish()= 0;

        NewmarkBeta<nbrdof_, nbrin_, nbrout_, real_type_> *system;
    };
    
    static constexpr int nbrdof= nbrdof_;
    static constexpr int nbrin= nbrin_;
    static constexpr int nbrout= nbrout_;

    typedef real_type_ real_type;
    typedef Eigen::Matrix<real_type, nbrdof_, 1> VecQ;
    typedef Eigen::Matrix<real_type, 2*nbrdof_, 1> VecX;
    typedef Eigen::Matrix<real_type, nbrin_, 1> VecI;
    typedef Eigen::Matrix<real_type, nbrout_, 1> VecO; 
    typedef Eigen::Matrix<real_type, nbrdof_, nbrdof_> MatQ;
    typedef Eigen::Matrix<real_type, nbrdof_, nbrin_> MatQI;
    typedef Eigen::Matrix<real_type, nbrout_, nbrdof_> MatOQ;
    typedef Eigen::Matrix<real_type, nbrout_, 2*nbrdof_+nbrin_> MatCD;
    typedef Eigen::Matrix<real_type, 2*nbrdof_, 2*nbrdof_> Mat2Q;
    typedef Eigen::Matrix<real_type, 2*nbrdof_, 2*nbrdof_ + nbrin_> MatS;
    
    NewmarkBeta(const std::string &aname= "anonymous_mbs", const std::string &adesc= "No description");
    virtual ~NewmarkBeta() {}

    void writeStateVariablesHeader(std::ostream &OutFile);
    void writeStateVariables(std::ostream &OutFile);
    VecQ computeResidualsInt();
    virtual VecQ computeResiduals()= 0;
    virtual void calcJacobian(real_type alphaM, real_type alphaC, real_type alphaK);
    virtual void calcB();
    virtual void calcCDF();
    virtual void calcOut()= 0;
    void staticEquilibrium();
    void staticEquilibriumWithLin();
    int newmarkOneStep(real_type h, bool hmodified= true);
    void newmarkInterval(real_type tfinal, real_type &h, real_type hmax);
    void newmarkSensitivities(real_type ts);
    void newmarkIntervalWithSens(real_type ts) { return newmarkIntervalWithSens(ts, ts); };
    void newmarkIntervalWithSens(real_type ts, real_type h);
    void newmarkIntegration(real_type tfinal, real_type hsave, real_type hmax, AbstractIntegratorVisitor *visitor= nullptr);
    void setOptionsFromFile(const std::string &fileName);

    std::string name;
    std::string description;

    std::array<std::string, nbrdof_> state_name;
    std::array<std::string, nbrin_> in_name;
    
    MatQ M;
    MatQ C;
    MatQ K;
    MatQI B;
    MatCD CD;
    MatOQ F;
    MatQ Jacobian;
    Eigen::FullPivLU<MatQ> LU;
    MatS S;
    
    VecX x;
    Eigen::Map<VecQ> q, qd;
    VecQ qdd;
    VecQ f;
    
    std::array<bool, nbrdof_> doflocked;
    VecI u;
    VecO y;
    
    real_type t= 0.0;
    real_type jac_fd_tol= 1e-2;
    real_type AbsTol= 1E-6;
    real_type RelTol= 1E-6;
    real_type StepTol= 1E-8;
    real_type hminmin= 1E-8;
    int jac_recalc_step= 4;
    int max_steps= 10;
    real_type Beta= 0.25;
    real_type Gamma=0.5;
    
    int n_back_steps;
    int n_steps;
    int n_sub_steps;
    real_type errq;
};

template <int nbrdof_, int nbrin_, int nbrout_, class real_type_>
NewmarkBeta<nbrdof_, nbrin_, nbrout_, real_type_>::NewmarkBeta(const std::string &aname, const std::string &adesc):
        name(aname),
        description(adesc),
        q(x.data()),
        qd(x.data()+nbrdof_)
    {
        x.setZero();
        qdd.setZero();
        doflocked.fill(false);
        u.setZero();
        for(int i= 0; i<nbrdof_; i++) {
            std::stringstream ss;
            ss << "q_" << i;
            state_name[i]= ss.str();
        }
        for(int i= 0; i<nbrin_; i++) {
            std::stringstream ss;
            ss << "u_" << i;
            in_name[i]= ss.str();
        }
    }

template <int nbrdof_, int nbrin_, int nbrout_, class real_type_>
void NewmarkBeta<nbrdof_, nbrin_, nbrout_, real_type_>::writeStateVariablesHeader(std::ostream &OutFile) {
    OutFile << " time ";

    for(int idof= 0; idof < nbrdof_; idof++)
        OutFile << " q" << idof << " qd" << idof << " qdd" << idof;

    for(int iinput= 0; iinput < nbrin_; iinput++)
        OutFile << " u" << iinput;

    OutFile << std::endl;
}

template <int nbrdof_, int nbrin_, int nbrout_, class real_type_>
void NewmarkBeta<nbrdof_, nbrin_, nbrout_, real_type_>::writeStateVariables(std::ostream &OutFile) {
    OutFile << " " << t;

    for(int idof= 0;idof < nbrdof_; idof++)
        OutFile << " " << q[idof] << " " << qd[idof] << " " << qdd[idof];

    for(int iinput= 0; iinput < nbrin_; iinput++)
        OutFile << " " << u[iinput];

    OutFile << std::endl;
}

template <int nbrdof_, int nbrin_, int nbrout_, class real_type_>
typename NewmarkBeta<nbrdof_, nbrin_, nbrout_, real_type_>::VecQ NewmarkBeta<nbrdof_, nbrin_, nbrout_, real_type_>::computeResidualsInt() {
    VecQ f= computeResiduals();
    for(int idof= 0; idof < nbrdof_; idof++)
        if(doflocked[idof]) f[idof]= 0.0;

    return f;
}

template <int nbrdof_, int nbrin_, int nbrout_, class real_type_>
void NewmarkBeta<nbrdof_, nbrin_, nbrout_, real_type_>::calcJacobian(real_type alphaM, real_type alphaC, real_type alphaK) {
    VecQ foff;
    f= computeResidualsInt();

    for(int jddl= 0; jddl < nbrdof_; jddl++) {
        real_type q_= q[jddl];
        real_type qd_= qd[jddl];
        real_type qdd_= qdd[jddl];

        q[jddl]+= alphaK*jac_fd_tol;
        qd[jddl]+= alphaC*jac_fd_tol;
        qdd[jddl]+= alphaM*jac_fd_tol;

        if(!doflocked[jddl]) {
            foff= computeResidualsInt();
            Jacobian.col(jddl)= (foff - f)/jac_fd_tol;
        } else {
            Jacobian.col(jddl).setZero();
            Jacobian(jddl, jddl)= 1.0;
        }

        q[jddl]= q_;
        qd[jddl]= qd_;
        qdd[jddl]= qdd_;
    }
}

template <int nbrdof_, int nbrin_, int nbrout_, class real_type_>
void NewmarkBeta<nbrdof_, nbrin_, nbrout_, real_type_>::calcB() {
    VecQ foff;
    f= computeResidualsInt();
    
    for(int jddl= 0; jddl < nbrin_; jddl++) {
        real_type u_= u[jddl];
        
        u[jddl]+= jac_fd_tol;
        
        foff= computeResidualsInt();
        B.col(jddl)= (foff - f)/jac_fd_tol;
        
        u[jddl]= u_;
    }
}

template <int nbrdof_, int nbrin_, int nbrout_, class real_type_>
void NewmarkBeta<nbrdof_, nbrin_, nbrout_, real_type_>::calcCDF() {
    calcOut();
    VecO y_base= y;
    
    for(int jddl= 0; jddl < nbrdof_; jddl++) {
        real_type q_= q[jddl];
        
        q[jddl]+= jac_fd_tol;
        
        calcOut();
        CD.col(jddl)= (y - y_base)/jac_fd_tol;
        
        q[jddl]= q_;
    }
    for(int jddl= 0; jddl < nbrdof_; jddl++) {
        real_type qd_= qd[jddl];
        
        qd[jddl]+= jac_fd_tol;
        
        calcOut();
        CD.col(jddl+nbrdof_)= (y - y_base)/jac_fd_tol;
        
        qd[jddl]= qd_;
    }
    for(int jddl= 0; jddl < nbrin_; jddl++) {
        real_type u_= u[jddl];
        
        u[jddl]+= jac_fd_tol;
        
        calcOut();
        CD.col(jddl+2*nbrdof_)= (y - y_base)/jac_fd_tol;
        
        u[jddl]= u_;
    }
    
    
    for(int jddl= 0; jddl < nbrdof_; jddl++) {
        real_type qdd_= qdd[jddl];
        
        qdd[jddl]+= jac_fd_tol;
        
        calcOut();
        F.col(jddl)= (y - y_base)/jac_fd_tol;
        
        qdd[jddl]= qdd_;
    }
    y= y_base;
}

template <int nbrdof_, int nbrin_, int nbrout_, class real_type_>
void NewmarkBeta<nbrdof_, nbrin_, nbrout_, real_type_>::staticEquilibrium() {
    real_type err;
    
    if(!is_finite(q) || !is_finite(qd) || !is_finite(qdd) || !is_finite(u))
        throw NewmarkBetaException("Static equilibrium calculation: initial values of q, qd, qdd, u are not finite.");

    int nstep= 0;
    do {
        nstep++;
        if((nstep%10)==1) {
            calcJacobian(100.0, 0.0, 1.0);
            LU.compute(Jacobian);
            if(!LU.isInvertible())
                throw NewmarkBetaException("Static equilibrium calculation: Jacobian is not invertible.");
        } else
            f= computeResidualsInt();

        VecQ qdd_corr= LU.solve(f);

        if(!is_finite(qdd_corr))
            throw NewmarkBetaException("Static equilibrium calculation: Newton step is not finite.");
            
        q-= qdd_corr;
        err= qdd_corr.norm();
    }
    while((nstep<1000) && (err > (1E-8 * sqrt(1.0*nbrdof_))));

    if(err > (1E-8 * sqrt(1.0*nbrdof_)))
        throw NewmarkBetaException("Static equilibrium calculation: tolerance not reached after 1000 iterations.");
    
}

template <int nbrdof_, int nbrin_, int nbrout_, class real_type_>
void NewmarkBeta<nbrdof_, nbrin_, nbrout_, real_type_>::staticEquilibriumWithLin() {
    staticEquilibrium();
    calcJacobian(0.0, 0.0, 0.0);
    calcB();
    
    // return continuous time linearization
    S.block(0, 0, nbrdof_, nbrdof_)= M;
    S.block(0, nbrdof_, nbrdof_, nbrdof_)= MatQ::Identity();
    S.block(nbrdof_, 0, nbrdof_, nbrdof_)= -C;
    S.block(nbrdof_, nbrdof_, nbrdof_, nbrdof_)= -K;
    
    S.block(0, 2*nbrdof_, nbrdof_, nbrin_).setZero();
    S.block(nbrdof_, 2*nbrdof_, nbrdof_, nbrin_)= -B;
}    
    
template <int nbrdof_, int nbrin_, int nbrout_, class real_type_>
int NewmarkBeta<nbrdof_, nbrin_, nbrout_, real_type_>::newmarkOneStep(real_type h, bool hmodified) {
    VecQ qdd_sto= qdd;
    t+= h;
    q+= h*qd + 0.5*h*h*qdd;
    qd+= h*qdd;

    int istepjac= 0;
    if(hmodified && jac_recalc_step>1) istepjac= 1;

    int nstep=0;
    real_type err;
    do {
        nstep++;
        n_sub_steps++;
        if((nstep%jac_recalc_step)==istepjac) {
            calcJacobian(1.0, Gamma*h, Beta*h*h);
            LU.compute(Jacobian);
            if(!LU.isInvertible())
                throw NewmarkBetaException("Step calculation: Jacobian is not invertible.");
        } else
            f= computeResidualsInt();

        VecQ qdd_corr= LU.solve(f);
        
        if(!is_finite(qdd_corr))
            throw NewmarkBetaException("Step calculation: Newton step is not finite.");
        
        err= qdd_corr.norm() / (sqrt(1.0*nbrdof_) * (1.0 + qdd.norm()));

        qdd-= qdd_corr;
        qdd_corr*= Gamma*h;
        qd-= qdd_corr;
        qdd_corr*= Beta*h/Gamma;
        q-= qdd_corr;
    } while((nstep<=max_steps) && (err>StepTol));

    if(nstep>max_steps) return 1;

    qdd_sto-= qdd;
    errq= 0.0;
    for(int iddl=0; iddl<nbrdof_; iddl++) {
        real_type errtmp, absqiddl;
        errtmp= h*h*fabs(qdd_sto[iddl]) / 12.0;
        absqiddl= fabs(q[iddl]);
        absqiddl= std::max(absqiddl, h*fabs(qd[iddl]));
        absqiddl= std::max(absqiddl, h*h*fabs(qdd[iddl]));
        if((absqiddl*RelTol)>AbsTol) errtmp= errtmp*AbsTol/(absqiddl*RelTol);
        errq+= errtmp*errtmp;
    }
    errq= sqrt(errq/nbrdof_);

    if (!std::isfinite(errq)) return 2;
    return 0;
}

template <int nbrdof_, int nbrin_, int nbrout_, class real_type_>
void NewmarkBeta<nbrdof_, nbrin_, nbrout_, real_type_>::newmarkInterval(real_type tfinal, real_type &h, real_type hmax) {
    bool hchanged= true;
    n_steps= 0;
    n_sub_steps= 0;
    n_back_steps= 0;

    if(!is_finite(q) || !is_finite(qd) || !is_finite(qdd) || !is_finite(u))
        throw NewmarkBetaException("Interval calculation: initial values of q, qd, qdd, u are not finite.");
    
    if (h>hmax) h= hmax;
    while(t < tfinal) {
        if((t+1.4*h) >= tfinal) {
            h= tfinal-t;
            hchanged= true;
        }
        real_type timesto= t;

        VecQ q_sto= q;
        VecQ qd_sto= qd;
        VecQ qdd_sto= qdd;

        int code= newmarkOneStep(h, hchanged);
        n_steps++;
        hchanged= false;

        if((code) || (errq>AbsTol))    {
            if((code==1) || (code==2) || !std::isfinite(errq))
                h*= 0.25;
            else
                h*= sqrt((0.21*AbsTol + 0.04*errq) / errq);
            hchanged= true;
            
            t= timesto;
            q= q_sto;
            qd= qd_sto;
            qdd= qdd_sto;
            
            if(h<hminmin)
                throw NewmarkBetaException("Interval calculation: dynamic step size too small.");
            
            n_back_steps++;
        } else {
            if((errq < (0.1*AbsTol)) && (h<hmax)) {
                h*= sqrt(AbsTol/(2.1*errq + 0.04*AbsTol));
                if(h>hmax) h= hmax;
                hchanged= true;
            }
        }
    }
}

template <int nbrdof_, int nbrin_, int nbrout_, class real_type_>
void NewmarkBeta<nbrdof_, nbrin_, nbrout_, real_type_>::newmarkIntervalWithSens(real_type ts, real_type h) {
    newmarkInterval(t+ts, h, ts);
    newmarkSensitivities(ts);
}

template <int nbrdof_, int nbrin_, int nbrout_, class real_type>
void NewmarkBeta<nbrdof_, nbrin_, nbrout_, real_type>::newmarkSensitivities(real_type ts) {
    // calculate pfk1_pqdd, pfk1_pqd, pfk1_pq (M, C, K), does not work with base class
    // TODO: maybe reuse Jacobian for next step or from last
    calcJacobian(1.0, Gamma*ts, Beta*ts*ts);
    calcB();
    calcCDF();
    
    MatQ pfk1_pqdk= C + ts*K;
    
    // pfk1_pqddk1= Jacobian
    LU.compute(Jacobian);
    if(!LU.isInvertible())
        throw NewmarkBetaException("Sensitivities calculation: Jacobian is not invertible.");
    
    MatQ pqddk1_pqk= -LU.solve(K);
    MatQ pqddk1_pqdk= -LU.solve(pfk1_pqdk);
    MatQI pqddk1_puk= -LU.solve(B);
    
    // pqk1_pqk
    S.block(0, 0, nbrdof_, nbrdof_)= MatQ::Identity() + ts*ts*Beta*pqddk1_pqk;
    // pqk1_pqdk
    S.block(0, nbrdof_, nbrdof_, nbrdof_)= ts*MatQ::Identity() + ts*ts*Beta*pqddk1_pqdk;
    // pqdk1_pqk
    S.block(nbrdof_, 0, nbrdof_, nbrdof_)= ts*Gamma*pqddk1_pqk;
    // pqdk1_pqdk
    S.block(nbrdof_, nbrdof_, nbrdof_, nbrdof_)= MatQ::Identity() + ts*Gamma*pqddk1_pqdk;

    // pqk1_puk
    S.block(0, 2*nbrdof_, nbrdof_, nbrin_)= ts*ts*Beta*pqddk1_puk;
    // pqdk1_puk
    S.block(nbrdof_, 2*nbrdof_, nbrdof_, nbrin_)= ts*Gamma*pqddk1_puk;
    
    CD.block(0, 0, nbrout_, nbrdof_)+= F*pqddk1_pqk;
    CD.block(0, nbrdof_, nbrout_, nbrdof_)+= F*pqddk1_pqdk;
    CD.block(0, 2*nbrdof_, nbrout_, nbrin_)+= F*pqddk1_puk;
}

template <int nbrdof_, int nbrin_, int nbrout_, class real_type_>
void NewmarkBeta<nbrdof_, nbrin_, nbrout_, real_type_>::newmarkIntegration(real_type tfinal, real_type hsave, real_type hmax, AbstractIntegratorVisitor *visitor) {
    t= 0.0;
    real_type h= 0.0;

    if(!is_finite(q) || !is_finite(qd) || !is_finite(qdd) || !is_finite(u))
        throw NewmarkBetaException("Integrator: initial values of q, qd, qdd, u are not finite.");

    if(visitor) {
        visitor->setSystem(this);
        visitor->start();
    }

    newmarkOneStep(0.0);
    if(visitor) visitor->step();

    int ipas= 0;
    bool res= true;
    h= hmax;
    while(t<tfinal)    {
        ipas++;
        newmarkInterval(hsave*ipas, h, hmax);
        
        if(visitor) visitor->step();
    }

    if(visitor) visitor->finish();
}

template <int nbrdof_, int nbrin_, int nbrout_, class real_type_>
void NewmarkBeta<nbrdof_, nbrin_, nbrout_, real_type_>::setOptionsFromFile(const std::string &fileName) {
    std::ifstream infile(fileName);
    real_type value;
    
    if(infile.is_open()) {
        int i= 0;
        for(std::string line; std::getline(infile, line);) {
            ++i;
            std::istringstream iss(line);
            
            std::string optionName;
            iss >> optionName;
            
            if(optionName.empty() || (optionName[0]=='*' || optionName[0]=='#')) // comment or empty line
                continue;
            
            iss >> value;
            if(iss.fail()) {
                fprintf(stderr, "Could not read value for option \"%s\" in line %d.\n", optionName.c_str(), i);
            }                
            
            if(optionName=="AbsTol") {
                if(value<=0.0) std::runtime_error("Option \"" + optionName + "\" must be positive.");
                AbsTol= value;
            } else if(optionName=="RelTol") {
                if(value<=0.0) std::runtime_error("Option \"" + optionName + "\" must be positive.");
                RelTol= value;
            } else if(optionName=="StepTol") {
                if(value<=0.0) std::runtime_error("Option \"" + optionName + "\" must be positive.");
                StepTol= value;
            } else if(optionName=="HMin") {
                if(value<=0.0) std::runtime_error("Option \"" + optionName + "\" must be positive.");
                hminmin= value;
            } else if(optionName=="JacRecalc") {
                if(value<1.0) std::runtime_error("Option \"" + optionName + "\" must be greater equal 1.");
                jac_recalc_step= value;
            } else if(optionName=="MaxSteps") {
                if(value<1.0) std::runtime_error("Option \"" + optionName + "\" must be greater equal 1.");
                max_steps= value;
            } else if(optionName=="Beta") {
                if(value<=0.0 || value >=1.0) std::runtime_error("Option \"" + optionName + "\" must be in (0, 1).");
                Beta= value;
            } else if(optionName=="Gamma") {
                if(value<=0.0 || value >=1.0) std::runtime_error("Option \"" + optionName + "\" must be in (0, 1).");
                Gamma= value;
            } else {
                std::runtime_error("Unknown Option \"" + optionName + "\".");
            }
        }            
    } else
        throw std::runtime_error("Could not open options file \"" + fileName + "\". Valid options are AbsTol, RelTol, StepTol, HMin, JacRecalc, MaxSteps, Beta, Gamma");
}   


#endif /* NEWMARKBETA_HPP_ */
