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
        Jacobian(nbrdof_, nbrdof_),
        LU(),
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

void ODEOrder2::calcJacobian(double alphaM, double alphaC, double alphaK, double tol) {
    VecX foff;
    f= computeResidualsInt();

    for(int jddl= 0; jddl < nbrdof; jddl++) {
        double q_= q[jddl];
        double qd_= qd[jddl];
        double qdd_= qdd[jddl];

        q[jddl]+= alphaK*tol;
        qd[jddl]+= alphaC*tol;
        qdd[jddl]+= alphaM*tol;

        if(!doflocked[jddl]) {
            foff= computeResidualsInt();
            Jacobian.col(jddl)= (foff - f)/tol;
        } else
            Jacobian.col(jddl).setZero();

        q[jddl]= q_;
        qd[jddl]= qd_;
        qdd[jddl]= qdd_;
    }
}

bool ODEOrder2::staticEquilibrium() {
    double err;

    int nstep= 0;
    do {
        nstep++;
        if((nstep%10)==1) {
            calcJacobian(100.0, 0.0, 1.0, 1.0E-5);
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

int ODEOrder2::newmarkOneStep(double h, bool hmodified) {
    const double Beta= 0.25;
    const double Gamma=0.5;

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
            calcJacobian(1.0, Gamma*h, Beta*h*h, 1E-2);
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
