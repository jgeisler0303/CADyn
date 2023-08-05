/*
 * NewmarkBeta_GMRES.hpp
 *
 *  Created on: 21.12.2022
 *      Author: jgeisler
 */

#ifndef NEWMARKBETA_GMRES_HPP_
#define NEWMARKBETA_GMRES_HPP_

#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <sstream>
#include <vector>
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <Eigen/IterativeLinearSolvers>
#include <unsupported/Eigen/IterativeSolvers>


template <int nbrdof_, int nbrin_, int nbrout_, class real_type, class ext_type>
class NewmarkBeta_GMRES;

namespace Eigen {
namespace internal {
// MatrixReplacement looks-like a Matrix, so let's inherits its traits:
template<int nbrdof_, int nbrin_, int nbrout_, class real_type, class ext_type>
struct traits<NewmarkBeta_GMRES<nbrdof_, nbrin_, nbrout_, real_type, ext_type>> :  public Eigen::internal::traits<Eigen::SparseMatrix<real_type> >
{};
}
}


template <int nbrdof_, int nbrin_, int nbrout_, class real_type, class ext_type>
class NewmarkBeta_GMRES : public Eigen::EigenBase<NewmarkBeta_GMRES<nbrdof_, nbrin_, nbrout_, real_type, ext_type>> {

public:
    static const int nbrdof;
    static const int nbrin;
    static const int nbrout;

    // Required typedefs, constants, and method:
    typedef real_type Scalar;
    typedef real_type RealScalar;
    typedef int StorageIndex;
    enum {
        ColsAtCompileTime = nbrdof,
        MaxColsAtCompileTime = nbrdof,
        IsRowMajor = false
    };
    
    typename Eigen::EigenBase<NewmarkBeta_GMRES<nbrdof_, nbrin_, nbrout_, real_type, ext_type>>::Index rows() const { return nbrdof; }
    typename Eigen::EigenBase<NewmarkBeta_GMRES<nbrdof_, nbrin_, nbrout_, real_type, ext_type>>::Index cols() const { return nbrdof; }
    
    template<typename Rhs>
    Eigen::Product<NewmarkBeta_GMRES<nbrdof_, nbrin_, nbrout_, real_type, ext_type>,Rhs,Eigen::AliasFreeProduct> operator*(const Eigen::MatrixBase<Rhs>& x) const {
        return Eigen::Product<NewmarkBeta_GMRES<nbrdof_, nbrin_, nbrout_, real_type, ext_type>,Rhs,Eigen::AliasFreeProduct>(*this, x.derived());
    }
    
    typedef Eigen::Matrix<real_type, nbrdof_, 1> VecQ;
    typedef Eigen::Matrix<real_type, nbrin_, 1> VecI;
    typedef Eigen::Matrix<real_type, nbrout_, 1> VecO; 
    typedef Eigen::Matrix<real_type, nbrdof_, nbrdof_> MatQ;
    typedef Eigen::Matrix<real_type, nbrdof_, nbrin_> MatQI;
    typedef Eigen::Matrix<real_type, nbrout_, nbrdof_> MatOQ;
    typedef Eigen::Matrix<real_type, nbrout_, 2*nbrdof_+nbrin_> MatCD;
    typedef Eigen::Matrix<real_type, 2*nbrdof_, 2*nbrdof_> Mat2Q;
    typedef Eigen::Matrix<real_type, 2*nbrdof_, 2*nbrdof_ + nbrin_> MatS;
    
    NewmarkBeta_GMRES(const std::string &aname= "anonymous_mbs", const std::string &adesc= "No description");
//     virtual ~NewmarkBeta_GMRES() {}

    VecQ computeResidualsInt(ext_type &ext, const VecQ &qdd, const VecQ &qd, const VecQ &q) const;
    virtual VecQ computeResiduals(ext_type &ext, const VecQ &qdd, const VecQ &qd, const VecQ &q) const;
    virtual VecQ calcJacobianProd(const VecQ &b) const;
    virtual void calcOut();
    
    bool staticEquilibrium();
    int newmarkOneStep(real_type h);
    bool newmarkInterval(real_type tfinal, real_type &h, real_type hmax);
    void setOptionsFromFile(const std::string &fileName);

    std::string name;
    std::string description;

    std::array<std::string, nbrdof_> state_name;
    std::array<std::string, nbrin_> in_name;
    
    Eigen::GMRES< NewmarkBeta_GMRES<nbrdof_, nbrin_, nbrout_, real_type, ext_type>, Eigen::IdentityPreconditioner > gmres;
    
    VecQ q, qd, qdd;
    VecQ f;
    
    std::array<bool, nbrdof_> doflocked;
    VecI u;
    VecO y;
    ext_type ext;
    
    real_type t= 0.0;
    real_type jac_fd_tol= 1e-2;
    real_type AbsTol= 1E-6;
    real_type RelTol= 1E-6;
    real_type StepTol= 1E-8;
    real_type hminmin= 1E-8;
    int jac_recalc_step= 4;
    int max_steps= 10;
    int gmres_iterations= 1;
    
    real_type Beta= 0.25;
    real_type Gamma=0.5;
    
    real_type alphaM;
    real_type alphaC;
    real_type alphaK;
    
    int n_back_steps;
    int n_steps;
    int n_sub_steps;
    real_type errq;
};

template <int nbrdof_, int nbrin_, int nbrout_, class real_type, class ext_type> const int NewmarkBeta_GMRES<nbrdof_, nbrin_, nbrout_, real_type, ext_type>::nbrdof= nbrdof_;
template <int nbrdof_, int nbrin_, int nbrout_, class real_type, class ext_type> const int NewmarkBeta_GMRES<nbrdof_, nbrin_, nbrout_, real_type, ext_type>::nbrin= nbrin_;
template <int nbrdof_, int nbrin_, int nbrout_, class real_type, class ext_type> const int NewmarkBeta_GMRES<nbrdof_, nbrin_, nbrout_, real_type, ext_type>::nbrout= nbrout_;

template <int nbrdof_, int nbrin_, int nbrout_, class real_type, class ext_type>
NewmarkBeta_GMRES<nbrdof_, nbrin_, nbrout_, real_type, ext_type>::NewmarkBeta_GMRES(const std::string &aname, const std::string &adesc):
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
        
        gmres.set_restart(nbrdof+1);
    }

template <int nbrdof_, int nbrin_, int nbrout_, class real_type, class ext_type>
typename NewmarkBeta_GMRES<nbrdof_, nbrin_, nbrout_, real_type, ext_type>::VecQ NewmarkBeta_GMRES<nbrdof_, nbrin_, nbrout_, real_type, ext_type>::computeResiduals(ext_type &ext, const VecQ &qdd, const VecQ &qd, const VecQ &q) const {
    return VecQ::Zero();
}

template <int nbrdof_, int nbrin_, int nbrout_, class real_type, class ext_type>
typename NewmarkBeta_GMRES<nbrdof_, nbrin_, nbrout_, real_type, ext_type>::VecQ NewmarkBeta_GMRES<nbrdof_, nbrin_, nbrout_, real_type, ext_type>::computeResidualsInt(ext_type &ext, const VecQ &qdd, const VecQ &qd, const VecQ &q) const {
    VecQ f= computeResiduals(ext, qdd, qd, q);
    for(int idof= 0; idof < nbrdof_; idof++)
        if(doflocked[idof]) f[idof]= 0.0;

    return f;
}

template <int nbrdof_, int nbrin_, int nbrout_, class real_type, class ext_type>
typename NewmarkBeta_GMRES<nbrdof_, nbrin_, nbrout_, real_type, ext_type>::VecQ NewmarkBeta_GMRES<nbrdof_, nbrin_, nbrout_, real_type, ext_type>::calcJacobianProd(const VecQ &b) const {
    VecQ df;

    VecQ q_= q;
    VecQ qd_= qd;
    VecQ qdd_= qdd;
    ext_type ext_;

    q_+= alphaK*jac_fd_tol * b;
    qd_+= alphaC*jac_fd_tol * b;
    qdd_+= alphaM*jac_fd_tol * b;

    df= computeResidualsInt(ext_, qdd_, qd_, q_);
    df= (df - f)/jac_fd_tol;

    return df;
}

template <int nbrdof_, int nbrin_, int nbrout_, class real_type, class ext_type>
bool NewmarkBeta_GMRES<nbrdof_, nbrin_, nbrout_, real_type, ext_type>::staticEquilibrium() {
    real_type err;

    alphaM= 100.0;
    alphaC= 0.0;
    alphaK= 1.0;
    
    gmres.compute(*this);
    gmres.setMaxIterations(nbrdof);
    
    int nstep= 0;
    do {
        nstep++;
        
        f= computeResidualsInt(ext, qdd, qd, q);
        
        Eigen::VectorXd b(nbrdof), x;
        b= f;
        x= gmres.solve(f);
        VecQ qdd_corr= x;

        q-= qdd_corr;
        err= qdd_corr.norm();
    }
    while((nstep<1000) && (err > (1E-8 * sqrt(1.0*nbrdof_))));

    return(err < (1E-8 * sqrt(1.0*nbrdof_)));
}

template <int nbrdof_, int nbrin_, int nbrout_, class real_type, class ext_type>
int NewmarkBeta_GMRES<nbrdof_, nbrin_, nbrout_, real_type, ext_type>::newmarkOneStep(real_type h) {
    VecQ qdd_sto= qdd;
    t+= h;
    q+= h*qd + 0.5*h*h*qdd;
    qd+= h*qdd;

    alphaM= 1.0;
    alphaC= Gamma*h;
    alphaK= Beta*h*h;
    
    gmres.compute(*this);
    gmres.setMaxIterations(gmres_iterations);

    int nstep=0;
    real_type err;
    do {
        nstep++;
        n_sub_steps++;
        
        f= computeResidualsInt(ext, qdd, qd, q);

        Eigen::VectorXd b(nbrdof), x;
        b= f;
        x= gmres.solve(f);
        VecQ qdd_corr= x;
//         VecQ qdd_corr= gmres.solve(f);
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

template <int nbrdof_, int nbrin_, int nbrout_, class real_type, class ext_type>
bool NewmarkBeta_GMRES<nbrdof_, nbrin_, nbrout_, real_type, ext_type>::newmarkInterval(real_type tfinal, real_type &h, real_type hmax) {
    n_steps= 0;
    n_sub_steps= 0;
    n_back_steps= 0;

    if (h>hmax || h<=0.0) h= hmax;
    while(t < tfinal) {
        if((t+1.4*h) >= tfinal) {
            h= tfinal-t;
        }
        real_type timesto= t;

        VecQ q_sto= q;
        VecQ qd_sto= qd;
        VecQ qdd_sto= qdd;

        int code= newmarkOneStep(h);
        
        n_steps++;
     
        if((code) || (errq>AbsTol))    {
            if((code==1) || (code==2) || !std::isfinite(errq))
                h*= 0.25;
            else
                h*= sqrt((0.21*AbsTol + 0.04*errq) / errq);
            
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
            }
        }
    }
    return true;
}

template <int nbrdof_, int nbrin_, int nbrout_, class real_type, class ext_type>
void NewmarkBeta_GMRES<nbrdof_, nbrin_, nbrout_, real_type, ext_type>::setOptionsFromFile(const std::string &fileName) {
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
            } else if(optionName=="GMRES_Iterations") {
                if(value<1.0) std::runtime_error("Option \"" + optionName + "\" must be greater equal 1.");
                gmres_iterations= value;
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

template <int nbrdof_, int nbrin_, int nbrout_, class real_type, class ext_type>
void NewmarkBeta_GMRES<nbrdof_, nbrin_, nbrout_, real_type, ext_type>::calcOut() {
}


// Implementation of NewmarkBeta_GMRES<nbrdof_, nbrin_, nbrout_, real_type, ext_type> * Eigen::DenseVector though a specialization of internal::generic_product_impl:
namespace Eigen {
namespace internal {
    template<typename Rhs, int nbrdof_, int nbrin_, int nbrout_, class real_type, class ext_type>
    struct generic_product_impl<NewmarkBeta_GMRES<nbrdof_, nbrin_, nbrout_, real_type, ext_type>, Rhs, SparseShape, DenseShape, GemvProduct> // GEMV stands for matrix-vector
    : generic_product_impl_base<NewmarkBeta_GMRES<nbrdof_, nbrin_, nbrout_, real_type, ext_type>, Rhs, generic_product_impl<NewmarkBeta_GMRES<nbrdof_, nbrin_, nbrout_, real_type, ext_type>,Rhs> >
    {
            typedef typename Product<NewmarkBeta_GMRES<nbrdof_, nbrin_, nbrout_, real_type, ext_type>, Rhs>::Scalar Scalar;
        
            template<typename Dest>
            static void scaleAndAddTo(Dest& dst, const NewmarkBeta_GMRES<nbrdof_, nbrin_, nbrout_, real_type, ext_type>& lhs, const Rhs& rhs, const Scalar& alpha)
            {
                // This method should implement "dst += alpha * lhs * rhs" inplace,
                // however, for iterative solvers, alpha is always equal to 1, so let's not bother about it.
                assert(alpha==Scalar(1) && "scaling is not implemented");
                EIGEN_ONLY_USED_FOR_DEBUG(alpha);
        
                // Here we could simply call dst.noalias() += lhs.my_matrix() * rhs,
                // but let's do something fancier (and less efficient):
                dst += lhs.calcJacobianProd(rhs);
            }; 
    };
}
}


#endif /* NEWMARKBETA_GMRES_HPP_ */
