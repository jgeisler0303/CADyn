/*
 * NewmarkBeta.hpp
 *
 *  Created on: 10.03.2017
 *      Author: jgeisler
 */

#ifndef NEWMARKBETA_HPP_
#define NEWMARKBETA_HPP_

#include <iosfwd>
#include <cmath>
#include <algorithm>
#include <sstream>
#include <vector>
#include <Eigen/Dense>
#include <Eigen/Geometry>


template <int nbrdof_, int nbrin_, int nbrout_, class real_type>
class NewmarkBeta {
public:
    class AbstractIntegratorVisitor {
    public:
        AbstractIntegratorVisitor(): system(nullptr) {}
        virtual ~AbstractIntegratorVisitor() {}

        void setSystem(NewmarkBeta<nbrdof_, nbrin_, nbrout_, real_type> *system_) { system= system_; }
        virtual void start()= 0;
        virtual void step()= 0;
        virtual void finish()= 0;

        NewmarkBeta<nbrdof_, nbrin_, nbrout_, real_type> *system;
    };
    
    static constexpr int nbrdof= nbrdof_;
    static constexpr int nbrin= nbrin_;
    static constexpr int nbrout= nbrout_;

    typedef Eigen::Matrix<real_type, nbrdof_, 1> VecQ;
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
    bool staticEquilibrium();
    bool staticEquilibriumWithLin();
    int newmarkOneStep(real_type h, bool hmodified= true);
    bool newmarkInterval(real_type tfinal, real_type &h, real_type hmax);
    bool newmarkIntervalWithSens(real_type ts) { return newmarkIntervalWithSens(ts, ts); };
    bool newmarkIntervalWithSens(real_type ts, real_type h);
    bool newmarkIntegration(real_type tfinal, real_type hsave, real_type hmax, AbstractIntegratorVisitor *visitor= nullptr);
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
    
    VecQ q, qd, qdd;
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
    real_type errq;
};

template <int nbrdof_, int nbrin_, int nbrout_, class real_type>
NewmarkBeta<nbrdof_, nbrin_, nbrout_, real_type>::NewmarkBeta(const std::string &aname, const std::string &adesc):
        name(aname),
        description(adesc)
    {
        q.setZero();
        qd.setZero();
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

template <int nbrdof_, int nbrin_, int nbrout_, class real_type>
void NewmarkBeta<nbrdof_, nbrin_, nbrout_, real_type>::writeStateVariablesHeader(std::ostream &OutFile) {
    OutFile << " time ";

    for(int idof= 0; idof < nbrdof_; idof++)
        OutFile << " q" << idof << " qd" << idof << " qdd" << idof;

    for(int iinput= 0; iinput < nbrin_; iinput++)
        OutFile << " u" << iinput;

    OutFile << std::endl;
}

template <int nbrdof_, int nbrin_, int nbrout_, class real_type>
void NewmarkBeta<nbrdof_, nbrin_, nbrout_, real_type>::writeStateVariables(std::ostream &OutFile) {
    OutFile << " " << t;

    for(int idof= 0;idof < nbrdof_; idof++)
        OutFile << " " << q[idof] << " " << qd[idof] << " " << qdd[idof];

    for(int iinput= 0; iinput < nbrin_; iinput++)
        OutFile << " " << u[iinput];

    OutFile << std::endl;
}

template <int nbrdof_, int nbrin_, int nbrout_, class real_type>
typename NewmarkBeta<nbrdof_, nbrin_, nbrout_, real_type>::VecQ NewmarkBeta<nbrdof_, nbrin_, nbrout_, real_type>::computeResidualsInt() {
    VecQ f= computeResiduals();
    for(int idof= 0; idof < nbrdof_; idof++)
        if(doflocked[idof]) f[idof]= 0.0;

    return f;
}

template <int nbrdof_, int nbrin_, int nbrout_, class real_type>
void NewmarkBeta<nbrdof_, nbrin_, nbrout_, real_type>::calcJacobian(real_type alphaM, real_type alphaC, real_type alphaK) {
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
        } else
            Jacobian.col(jddl).setZero();

        q[jddl]= q_;
        qd[jddl]= qd_;
        qdd[jddl]= qdd_;
    }
}

template <int nbrdof_, int nbrin_, int nbrout_, class real_type>
void NewmarkBeta<nbrdof_, nbrin_, nbrout_, real_type>::calcB() {
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

template <int nbrdof_, int nbrin_, int nbrout_, class real_type>
void NewmarkBeta<nbrdof_, nbrin_, nbrout_, real_type>::calcCDF() {
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

template <int nbrdof_, int nbrin_, int nbrout_, class real_type>
bool NewmarkBeta<nbrdof_, nbrin_, nbrout_, real_type>::staticEquilibrium() {
    real_type err;

    int nstep= 0;
    do {
        nstep++;
        if((nstep%10)==1) {
            calcJacobian(100.0, 0.0, 1.0);
            LU.compute(Jacobian);
        } else
            f= computeResidualsInt();

        VecQ qdd_corr= LU.solve(f);
        q-= qdd_corr;
        err= qdd_corr.norm();
    }
    while((nstep<1000) && (err > (1E-8 * sqrt(1.0*nbrdof_))));

    return(nstep<=1000);
}

template <int nbrdof_, int nbrin_, int nbrout_, class real_type>
bool NewmarkBeta<nbrdof_, nbrin_, nbrout_, real_type>::staticEquilibriumWithLin() {
    bool res= staticEquilibrium();
    calcJacobian(0.0, 0.0, 0.0);
    calcB();
    
    // MatX Pxn1_Pxn
    S.block(0, 0, nbrdof_, nbrdof_)= M;
    // MatX Pxn1_Pdxn
    S.block(0, nbrdof_, nbrdof_, nbrdof_)= MatQ::Identity();
    // MatX Pdxn1_Pxn
    S.block(nbrdof_, 0, nbrdof_, nbrdof_)= -C;
    // MatX Pdxn1_Pdxn
    S.block(nbrdof_, nbrdof_, nbrdof_, nbrdof_)= -K;
    
    S.block(0, 2*nbrdof_, nbrdof_, nbrin_).setZero();
    S.block(nbrdof_, 2*nbrdof_, nbrdof_, nbrin_)= -B;
    
    return res;
}    
    
template <int nbrdof_, int nbrin_, int nbrout_, class real_type>
int NewmarkBeta<nbrdof_, nbrin_, nbrout_, real_type>::newmarkOneStep(real_type h, bool hmodified) {
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
        if((nstep%jac_recalc_step)==istepjac) {
            calcJacobian(1.0, Gamma*h, Beta*h*h);
            LU.compute(Jacobian);
        } else
            f= computeResidualsInt();

        VecQ qdd_corr= LU.solve(f);
        err= qdd_corr.norm() / (sqrt(1.0*nbrdof_) * (1.0 + qdd.norm()));

        qdd-= qdd_corr;
        qdd_corr*= Gamma*h;
        qd-= qdd_corr;
        qdd_corr*= Beta*h/Gamma;
        q-= qdd_corr;
    } while((nstep<max_steps) && (err>StepTol));

    if(nstep==max_steps) return 1;

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

template <int nbrdof_, int nbrin_, int nbrout_, class real_type>
bool NewmarkBeta<nbrdof_, nbrin_, nbrout_, real_type>::newmarkInterval(real_type tfinal, real_type &h, real_type hmax) {
    bool hchanged= true;
    n_steps= 0;
    n_back_steps= 0;

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
                return false;
            
            n_back_steps++;
        } else {
            if((errq < (0.1*AbsTol)) && (h<hmax)) {
                h*= sqrt(AbsTol/(2.1*errq + 0.04*AbsTol));
                if(h>hmax) h= hmax;
                hchanged= true;
            }
        }
    }
    return true;
}

template <int nbrdof_, int nbrin_, int nbrout_, class real_type>
bool NewmarkBeta<nbrdof_, nbrin_, nbrout_, real_type>::newmarkIntervalWithSens(real_type ts, real_type h) {
    bool hchanged= true;
    bool first_run= true;
    n_steps= 0;
    n_back_steps= 0;
    real_type tfinal= t+ts;
    
    Mat2Q PG_Pddxn_ddxn1();
    MatS PG_Pxn_dxn_u();
    
    newmarkOneStep(0.0);
    
    while(t < tfinal) {
        if((t+1.4*h) >= tfinal) {
            h= tfinal-t;
            hchanged= true;
        }
        real_type timesto= t;

        VecQ q_sto= q;
        VecQ qd_sto= qd;
        VecQ qdd_sto= qdd;
        
        if(first_run) {
            calcJacobian(1.0, Gamma*h, Beta*h*h);
            calcB();
            calcCDF();
            PG_Pddxn_ddxn1.block(0, 0, nbrdof_, nbrdof_)= M;
            PG_Pxn_dxn_u.block(0, nbrdof_, nbrdof_, nbrdof_)= C;
            PG_Pxn_dxn_u.block(0, 0, nbrdof_, nbrdof_)= K;
            PG_Pxn_dxn_u.block(0, 2*nbrdof_, nbrdof_, nbrin_)= B;
            
            LU.compute(Jacobian);
            hchanged= false;
            first_run= false;
        }
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
                return false;
            
            n_back_steps++;
        } else {
            if((errq < (0.1*AbsTol)) && (h<ts)) {
                h*= sqrt(AbsTol/(2.1*errq + 0.04*AbsTol));
                if(h>ts) h= ts;
                hchanged= true;
            }
        }
    }
    
    calcJacobian(0.0, 0.0, 0.0);
    calcB();
    
    PG_Pxn_dxn_u.block(nbrdof_, 0, nbrdof_, nbrdof_)= K;
    PG_Pxn_dxn_u.block(nbrdof_, nbrdof_, nbrdof_, nbrdof_)= C + ts*K;
    PG_Pxn_dxn_u.block(nbrdof_, 2*nbrdof_, nbrdof_, nbrin_)= B;
    
    PG_Pddxn_ddxn1.block(nbrdof_, 0, nbrdof_, nbrdof_)= ts*(1.0-Gamma)*C + ts*ts*(0.5-Beta)*K;
    PG_Pddxn_ddxn1.block(0, nbrdof_, nbrdof_, nbrdof_)= MatQ::Zero();
    PG_Pddxn_ddxn1.block(nbrdof_, nbrdof_, nbrdof_, nbrdof_)= M + ts*Gamma*C + ts*ts*Beta*K;
    
    Eigen::ColPivHouseholderQR<Mat2Q> invPG_Pddxn_ddxn1(PG_Pddxn_ddxn1);
    // alternative: compute inverse via schure complement
    
    MatS Pddx_Pxn_dxn= -invPG_Pddxn_ddxn1.solve(PG_Pxn_dxn_u);
    
    #define Pddxn_Pxn   Pddx_Pxn_dxn.block(0, 0, nbrdof_, nbrdof_)
    #define Pddxn1_Pxn  Pddx_Pxn_dxn.block(nbrdof_, 0, nbrdof_, nbrdof_)
    #define Pddxn_Pdxn  Pddx_Pxn_dxn.block(0, nbrdof_, nbrdof_, nbrdof_)
    #define Pddxn1_Pdxn Pddx_Pxn_dxn.block(nbrdof_, nbrdof_, nbrdof_, nbrdof_)
    #define Pddxn_Pu  Pddx_Pxn_dxn.block(0, 2*nbrdof_, nbrdof_, nbrin_)
    #define Pddxn1_Pu Pddx_Pxn_dxn.block(nbrdof_, 2*nbrdof_, nbrdof_, nbrin_)
    
    // MatX Pxn1_Pxn
    S.block(0, 0, nbrdof_, nbrdof_)= MatQ::Identity() + ts*ts*((0.5-Beta)*Pddxn_Pxn + Beta*Pddxn1_Pxn);
    // MatX Pxn1_Pdxn
    S.block(0, nbrdof_, nbrdof_, nbrdof_)= ts*MatQ::Identity() + ts*ts*((0.5-Beta)*Pddxn_Pdxn + Beta*Pddxn1_Pdxn);
    // MatX Pdxn1_Pxn
    S.block(nbrdof_, 0, nbrdof_, nbrdof_)= ts*((1.0-Gamma)*Pddxn_Pxn + Gamma*Pddxn1_Pxn);
    // MatX Pdxn1_Pdxn
    S.block(nbrdof_, nbrdof_, nbrdof_, nbrdof_)= MatQ::Identity() + ts*((1.0-Gamma)*Pddxn_Pdxn + Gamma*Pddxn1_Pdxn);

    S.block(0, 2*nbrdof_, nbrdof_, nbrin_)= ts*ts*((0.5-Beta)*Pddxn_Pu + Beta*Pddxn1_Pu);
    S.block(nbrdof_, 2*nbrdof_, nbrdof_, nbrin_)= ts*((1.0-Gamma)*Pddxn_Pu + Gamma*Pddxn1_Pu);
    
    CD.block(0, 0, nbrout_, nbrdof_)+= F*Pddxn1_Pxn; // TODO: really ddxn_1_ ?
    CD.block(0, nbrdof_, nbrout_, nbrdof_)+= F*Pddxn1_Pdxn;
    CD.block(0, 2*nbrdof_, nbrout_, nbrin_)+= F*Pddxn1_Pu;
    
    return true;
}

template <int nbrdof_, int nbrin_, int nbrout_, class real_type>
bool NewmarkBeta<nbrdof_, nbrin_, nbrout_, real_type>::newmarkIntegration(real_type tfinal, real_type hsave, real_type hmax, AbstractIntegratorVisitor *visitor) {
    t= 0.0;
    real_type h= 0.0;

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
        if(!newmarkInterval(hsave*ipas, h, hmax)) {
            res= false;
            break;
        }
        if(visitor) visitor->step();
    }

    if(visitor) visitor->finish();
    return res;
}

template <int nbrdof_, int nbrin_, int nbrout_, class real_type>
void NewmarkBeta<nbrdof_, nbrin_, nbrout_, real_type>::setOptionsFromFile(const std::string &fileName) {
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
