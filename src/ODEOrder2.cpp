/*
 * ODEOrder2.cpp
 *
 *  Created on: 10.03.2017
 *      Author: jgeisler
 */

#include "ODEOrder2.hpp"
#include <cmath>
#include <algorithm>
#include <sstream>

ODEOrder2::ODEOrder2(int nbrdof_, int nbrin_, const string &aname, const string &adesc):
        name(aname),
        description(adesc),
        nbrdof(nbrdof_),
        nbrin(nbrin_),
        state_name(nbrdof_),
        in_name(nbrin_),
        M(nbrdof_, nbrdof_),
        C(nbrdof_, nbrdof_),
        K(nbrdof_, nbrdof_),
        B(nbrdof_, nbrin_),    
        Jacobian(nbrdof_, nbrdof_),
        LU(),
        S(2*nbrdof_, 2*nbrdof_ + nbrin_),
        q(nbrdof_),
        qd(nbrdof_),
        qdd(nbrdof_),
        f(nbrdof_),
        doflocked(nbrdof_, false),
        u(nbrin_),
        t(0.0)
    {
        q.setZero();
        qd.setZero();
        qdd.setZero();
        u.setZero();
        for(int i= 0; i<nbrdof; i++) {
            std::stringstream ss;
            ss << "q_" << i;
            state_name[i]= ss.str();
        }
        for(int i= 0; i<nbrin; i++) {
            std::stringstream ss;
            ss << "u_" << i;
            in_name[i]= ss.str();
        }
    }

void ODEOrder2::writeStateVariablesHeader(std::ostream &OutFile) {
    OutFile << " time ";

    for(int idof= 0; idof < nbrdof; idof++)
        OutFile << " q" << idof << " qd" << idof << " qdd" << idof;

    for(int iinput= 0; iinput < nbrin; iinput++)
        OutFile << " u" << iinput;

    OutFile << std::endl;
}

void ODEOrder2::writeStateVariables(std::ostream &OutFile) {
    OutFile << " " << t;

    for(int idof= 0;idof < nbrdof; idof++)
        OutFile << " " << q[idof] << " " << qd[idof] << " " << qdd[idof];

    for(int iinput= 0; iinput < nbrin; iinput++)
        OutFile << " " << u[iinput];

    OutFile << std::endl;
}

VecX ODEOrder2::computeResidualsInt() {
    VecX f= computeResiduals();
    for(int idof= 0; idof < nbrdof; idof++)
        if(doflocked[idof]) f[idof]= 0.0;

    return f;
}

void ODEOrder2::calcJacobian(double alphaM, double alphaC, double alphaK) {
    VecX foff;
    f= computeResidualsInt();

    for(int jddl= 0; jddl < nbrdof; jddl++) {
        double q_= q[jddl];
        double qd_= qd[jddl];
        double qdd_= qdd[jddl];

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

void ODEOrder2::calcB() {
    VecX foff;
    f= computeResidualsInt();
    
    for(int jddl= 0; jddl < nbrin; jddl++) {
        double u_= u[jddl];
        
        u[jddl]+= jac_fd_tol;
        
        foff= computeResidualsInt();
        B.col(jddl)= (foff - f)/jac_fd_tol;
        
        u[jddl]= u_;
    }
}

bool ODEOrder2::staticEquilibrium() {
    double err;

    int nstep= 0;
    do {
        nstep++;
        if((nstep%10)==1) {
            calcJacobian(100.0, 0.0, 1.0);
            LU.compute(Jacobian);
        } else
            f= computeResidualsInt();

        VecX qdd_corr= LU.solve(f);
        q-= qdd_corr;
        err= qdd_corr.norm();
    }
    while((nstep<1000) && (err > (1E-8 * sqrt(1.0*nbrdof))));

    return(nstep<=1000);
}

bool ODEOrder2::staticEquilibriumWithLin() {
    bool res= staticEquilibrium();
    calcJacobian(0.0, 0.0, 0.0);
    calcB();
    
    // MatX Pxn1_Pxn
    S.block(0, 0, nbrdof, nbrdof)= M;
    // MatX Pxn1_Pdxn
    S.block(0, nbrdof, nbrdof, nbrdof)= MatX::Identity(nbrdof, nbrdof);
    // MatX Pdxn1_Pxn
    S.block(nbrdof, 0, nbrdof, nbrdof)= -C;
    // MatX Pdxn1_Pdxn
    S.block(nbrdof, nbrdof, nbrdof, nbrdof)= -K;
    
    S.block(0, 2*nbrdof, nbrdof, nbrin).setZero();
    S.block(nbrdof, 2*nbrdof, nbrdof, nbrin)= -B;
    
    return res;
}    
    
int ODEOrder2::newmarkOneStep(double h, bool hmodified) {
    VecX qdd_sto= qdd;
    t+= h;
    q+= h*qd + 0.5*h*h*qdd;
    qd+= h*qdd;

    int istepjac= 0;
    if(hmodified) istepjac= 1;

    int nstep=0;
    double err;
    do {
        nstep++;
        if((nstep%jac_recalc_step)==istepjac) {
            calcJacobian(1.0, Gamma*h, Beta*h*h);
            LU.compute(Jacobian);
        } else
            f= computeResidualsInt();

        VecX qdd_corr= LU.solve(f);
        err= qdd_corr.norm() / (sqrt(1.0*nbrdof) * (1.0 + qdd.norm()));

        qdd-= qdd_corr;
        qdd_corr*= Gamma*h;
        qd-= qdd_corr;
        qdd_corr*= Beta*h/Gamma;
        q-= qdd_corr;
    } while((nstep<max_steps) && (err>StepTol));

    if(nstep==max_steps) return 1;

    qdd_sto-= qdd;
    errq= 0.0;
    for(int iddl=0; iddl<nbrdof; iddl++) {
        double errtmp, absqiddl;
        errtmp= h*h*fabs(qdd_sto[iddl]) / 12.0;
        absqiddl= fabs(q[iddl]);
        absqiddl= std::max(absqiddl, h*fabs(qd[iddl]));
        absqiddl= std::max(absqiddl, h*h*fabs(qdd[iddl]));
        if((absqiddl*RelTol)>AbsTol) errtmp= errtmp*AbsTol/(absqiddl*RelTol);
        errq+= errtmp*errtmp;
    }
    errq= sqrt(errq/nbrdof);

    if (!std::isfinite(errq)) return 2;
    return 0;
}

bool ODEOrder2::newmarkInterval(double tfinal, double &h, double hmax) {
    bool hchanged= true;
    n_steps= 0;
    n_back_steps= 0;

    if (h>hmax) h= hmax;
    while(t < tfinal) {
        if((t+1.4*h) >= tfinal) {
            h= tfinal-t;
            hchanged= true;
        }
        double timesto= t;

        VecX q_sto= q;
        VecX qd_sto= qd;
        VecX qdd_sto= qdd;

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

bool ODEOrder2::newmarkIntervalWithSens(double ts, double h) {
    bool hchanged= true;
    bool first_run= true;
    n_steps= 0;
    n_back_steps= 0;
    double tfinal= t+ts;
    
    MatX PG_Pddxn_ddxn1(2*nbrdof, 2*nbrdof);
    MatX PG_Pxn_dxn_u(2*nbrdof, 2*nbrdof + nbrin);
    
    newmarkOneStep(0.0);
    
    while(t < tfinal) {
        if((t+1.4*h) >= tfinal) {
            h= tfinal-t;
            hchanged= true;
        }
        double timesto= t;

        VecX q_sto= q;
        VecX qd_sto= qd;
        VecX qdd_sto= qdd;
        
        if(first_run) {
            calcJacobian(1.0, Gamma*h, Beta*h*h);
            calcB();
            PG_Pddxn_ddxn1.block(0, 0, nbrdof, nbrdof)= M;
            PG_Pxn_dxn_u.block(0, nbrdof, nbrdof, nbrdof)= C;
            PG_Pxn_dxn_u.block(0, 0, nbrdof, nbrdof)= K;
            PG_Pxn_dxn_u.block(0, 2*nbrdof, nbrdof, nbrin)= B;
            
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
    
    PG_Pxn_dxn_u.block(nbrdof, 0, nbrdof, nbrdof)= K;
    PG_Pxn_dxn_u.block(nbrdof, nbrdof, nbrdof, nbrdof)= C + ts*K;
    PG_Pxn_dxn_u.block(nbrdof, 2*nbrdof, nbrdof, nbrin)= B;
    
    PG_Pddxn_ddxn1.block(nbrdof, 0, nbrdof, nbrdof)= ts*(1.0-Gamma)*C + ts*ts*(0.5-Beta)*K;
    PG_Pddxn_ddxn1.block(0, nbrdof, nbrdof, nbrdof)= MatX::Zero(nbrdof, nbrdof);
    PG_Pddxn_ddxn1.block(nbrdof, nbrdof, nbrdof, nbrdof)= M + ts*Gamma*C + ts*ts*Beta*K;
    
    Eigen::ColPivHouseholderQR<MatX> invPG_Pddxn_ddxn1(PG_Pddxn_ddxn1);
    // alternative: compute inverse via schure complement
    
    MatX Pddx_Pxn_dxn= -invPG_Pddxn_ddxn1.solve(PG_Pxn_dxn_u);
    
    #define Pddxn_Pxn   Pddx_Pxn_dxn.block(0, 0, nbrdof, nbrdof)
    #define Pddxn1_Pxn  Pddx_Pxn_dxn.block(nbrdof, 0, nbrdof, nbrdof)
    #define Pddxn_Pdxn  Pddx_Pxn_dxn.block(0, nbrdof, nbrdof, nbrdof)
    #define Pddxn1_Pdxn Pddx_Pxn_dxn.block(nbrdof, nbrdof, nbrdof, nbrdof)
    #define Pddxn_Pu  Pddx_Pxn_dxn.block(0, 2*nbrdof, nbrdof, nbrin)
    #define Pddxn1_Pu Pddx_Pxn_dxn.block(nbrdof, 2*nbrdof, nbrdof, nbrin)
    
    // MatX Pxn1_Pxn
    S.block(0, 0, nbrdof, nbrdof)= MatX::Identity(nbrdof, nbrdof) + ts*ts*((0.5-Beta)*Pddxn_Pxn + Beta*Pddxn1_Pxn);
    // MatX Pxn1_Pdxn
    S.block(0, nbrdof, nbrdof, nbrdof)= ts*MatX::Identity(nbrdof, nbrdof) + ts*ts*((0.5-Beta)*Pddxn_Pdxn + Beta*Pddxn1_Pdxn);
    // MatX Pdxn1_Pxn
    S.block(nbrdof, 0, nbrdof, nbrdof)= ts*((1.0-Gamma)*Pddxn_Pxn + Gamma*Pddxn1_Pxn);
    // MatX Pdxn1_Pdxn
    S.block(nbrdof, nbrdof, nbrdof, nbrdof)= MatX::Identity(nbrdof, nbrdof) + ts*((1.0-Gamma)*Pddxn_Pdxn + Gamma*Pddxn1_Pdxn);

    S.block(0, 2*nbrdof, nbrdof, nbrin)= ts*ts*((0.5-Beta)*Pddxn_Pu + Beta*Pddxn1_Pu);
    S.block(nbrdof, 2*nbrdof, nbrdof, nbrin)= ts*((1.0-Gamma)*Pddxn_Pu + Gamma*Pddxn1_Pu);
    
    return true;
}

bool ODEOrder2::newmarkIntegration(double tfinal, double hsave, double hmax, AbstractIntegratorVisitor *visitor) {
    t= 0.0;
    double h= 0.0;

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
