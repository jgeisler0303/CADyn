#ifndef MULTIBODY_HPP_
#define MULTIBODY_HPP_

#include <string>
#include <array>
#include <Eigen/Dense>
#include <Eigen/Geometry>

template<typename scalar_type>
using Vec3= Eigen::Matrix<scalar_type, 3, 1>;

template<typename scalar_type>
using Mat3= Eigen::Matrix<scalar_type, 3, 3>;

template<typename scalar_type>
using Trans= Eigen::Transform<scalar_type, 3, Eigen::Affine>;

template<typename scalar_type, int nbrdof_>
using VecQ= Eigen::Matrix<scalar_type, nbrdof_, 1>;

template<typename scalar_type>
using Vec6= Eigen::Matrix<scalar_type, 6, 1>;

template<typename scalar_type, int nbredof_>
using vector_VecE= std::array<Vec3<scalar_type>, nbredof_>;

template<typename scalar_type, int nbredof_>
using  vector_MatE= std::array<Mat3<scalar_type>, nbredof_>;

template<typename scalar_type, int nbredof_>
using VecE= Eigen::Matrix<scalar_type, nbredof_, 1>;

template<typename scalar_type, int nbredof_>
using MatE= Eigen::Matrix<scalar_type, nbredof_, nbredof_>;

template<typename scalar_type, int nbredof_>
using MatE6= Eigen::Matrix<scalar_type, nbredof_, 6>;



template<typename scalar_type, typename real_type, int nbrdof_>
class Kinematics {
public:
    Kinematics() {
        vGpartial.setZero();
        omegapartial.setZero();
        Fext.setZero();
        Mext.setZero();
    }

    Trans<scalar_type> TG;
    Vec3<scalar_type> vG;
    Vec3<scalar_type> aG;
    Vec3<scalar_type> omega;
    Vec3<scalar_type> omegad;
    Vec3<scalar_type> R;
    Vec3<scalar_type> MG;
    Vec3<scalar_type> Fext;
    Vec3<scalar_type> Mext;

    Eigen::Matrix<scalar_type, 3, nbrdof_> vGpartial;
    Eigen::Matrix<scalar_type, 3, nbrdof_> omegapartial;
    
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

template<typename scalar_type, typename real_type, int nbrdof_>
class Body: public Kinematics<scalar_type, real_type, nbrdof_> {
public:
    using Kinematics<scalar_type, real_type, nbrdof_>::TG;
    using Kinematics<scalar_type, real_type, nbrdof_>::vG;
    using Kinematics<scalar_type, real_type, nbrdof_>::aG;
    using Kinematics<scalar_type, real_type, nbrdof_>::omega;
    using Kinematics<scalar_type, real_type, nbrdof_>::omegad;
    using Kinematics<scalar_type, real_type, nbrdof_>::R;
    using Kinematics<scalar_type, real_type, nbrdof_>::MG;
    using Kinematics<scalar_type, real_type, nbrdof_>::Fext;
    using Kinematics<scalar_type, real_type, nbrdof_>::Mext;

    using Kinematics<scalar_type, real_type, nbrdof_>::vGpartial;
    using Kinematics<scalar_type, real_type, nbrdof_>::omegapartial;
    
    Body(Kinematics<scalar_type, real_type, nbrdof_> *rel_base, const std::string &aname= "Body", const std::string &adesc= "No description"):
        name(aname),
        description(adesc),
        mass(0.0),
        PhiG(),
        rel_base(rel_base)
        {}

    std::string name;
    std::string description;

    void composeMotion();
    void computeGravity(const Vec3<real_type> &g);
    void computeForceBalance();
    void generalizedBalance(VecQ<scalar_type, nbrdof_> &f);

    real_type mass;
    Mat3<real_type> PhiG;

    Kinematics<scalar_type, real_type, nbrdof_> *rel_base;

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

template<typename scalar_type, typename real_type, int nbrdof_>
void Body<scalar_type, real_type, nbrdof_>::composeMotion() {
    Vec3<scalar_type> erel= rel_base->TG.linear() * TG.translation();
    Vec3<scalar_type> vrel= rel_base->TG.linear() * vG;
    Vec3<scalar_type> arel= rel_base->TG.linear() * aG;
    Vec3<scalar_type> wrel= rel_base->TG.linear() * omega;
    Vec3<scalar_type> wdrel=rel_base->TG.linear() * omegad;

    TG= rel_base->TG * TG;
    vG= rel_base->vG + rel_base->omega.cross(erel) + vrel;
    omega= rel_base->omega + wrel;
    aG= rel_base->aG + rel_base->omegad.cross(erel) + rel_base->omega.cross(rel_base->omega.cross(erel)) + 2.0*rel_base->omega.cross(vrel) + arel;
    omegad= rel_base->omegad + rel_base->omega.cross(wrel) + wdrel;

    for(int idof= 0; idof<vGpartial.cols(); idof++) {
        vrel= rel_base->TG.linear() * vGpartial.col(idof);
        wrel= rel_base->TG.linear() * omegapartial.col(idof);
        vGpartial.col(idof)= rel_base->vGpartial.col(idof) + rel_base->omegapartial.col(idof).cross(erel) + vrel;
        omegapartial.col(idof)= rel_base->omegapartial.col(idof) + wrel;
    }
}

template<typename scalar_type, typename real_type, int nbrdof_>
void Body<scalar_type, real_type, nbrdof_>::computeForceBalance() {
    R= Fext;
    R-=mass*aG;

    Mat3<scalar_type> PhiG_local= TG.linear()*PhiG*TG.linear().transpose();
    MG= Mext;
    MG-= PhiG_local * omegad;
    MG-= omega.cross(PhiG_local*omega);
}

template<typename scalar_type, typename real_type, int nbrdof_>
void Body<scalar_type, real_type, nbrdof_>::generalizedBalance(VecQ<scalar_type, nbrdof_> &f) {
    f-= vGpartial.transpose() * R;
    f-= omegapartial.transpose() * MG;
}

template<typename scalar_type, typename real_type, int nbrdof_>
void Body<scalar_type, real_type, nbrdof_>::computeGravity(const Vec3<real_type> &g) {
	R+= mass * g;
}


template <typename scalar_type, typename real_type, int nbrdof_, int nbredof_>
class ElasticBody: public Body<scalar_type, real_type, nbrdof_> {
public:
    using Kinematics<scalar_type, real_type, nbrdof_>::TG;
    using Kinematics<scalar_type, real_type, nbrdof_>::vG;
    using Kinematics<scalar_type, real_type, nbrdof_>::aG;
    using Kinematics<scalar_type, real_type, nbrdof_>::omega;
    using Kinematics<scalar_type, real_type, nbrdof_>::omegad;
    using Kinematics<scalar_type, real_type, nbrdof_>::R;
    using Kinematics<scalar_type, real_type, nbrdof_>::MG;
    using Kinematics<scalar_type, real_type, nbrdof_>::Fext;
    using Kinematics<scalar_type, real_type, nbrdof_>::Mext;

    using Kinematics<scalar_type, real_type, nbrdof_>::vGpartial;
    using Kinematics<scalar_type, real_type, nbrdof_>::omegapartial;

    using Body<scalar_type, real_type, nbrdof_>::mass;
    using Body<scalar_type, real_type, nbrdof_>::PhiG;

    using Body<scalar_type, real_type, nbrdof_>::rel_base;

    static constexpr int nbredof= nbredof_;

    ElasticBody(Kinematics<scalar_type, real_type, nbrdof_> *rel_base, std::array<int, nbredof_> edof_, const std::string &name_= "ElasticBody", const std::string &desc_= "No description"):
        Body<scalar_type, real_type, nbrdof_>(rel_base, name_, desc_),
        edof(edof_)
    {
        Fe_ext.setZero();
    }

    void computeGravity(const Vec3<real_type> &g);
    void computeForceBalance();
    void generalizedBalance(VecQ<scalar_type, nbrdof_> &f);
    void setEDOF(const VecQ<scalar_type, nbrdof_> &q, const VecQ<scalar_type, nbrdof_> &qd, const VecQ<scalar_type, nbrdof_> &qdd);

    std::array<int, nbredof_> edof;
    VecE<scalar_type, nbredof_> eq;
    VecE<scalar_type, nbredof_> eqd;
    VecE<scalar_type, nbredof_> eqdd;
    Vec3<real_type> md;
    vector_VecE<real_type, nbredof_> Ct;
    vector_VecE<real_type, nbredof_> Cr;
    MatE<real_type, nbredof_> Me;

    vector_MatE<real_type, nbredof_> Gr;
    vector_VecE<real_type, nbredof_> Ge;
    MatE6<real_type, nbredof_> Oe;

    VecE<scalar_type, nbredof_> Re;
    VecE<real_type, nbredof_> K;
    VecE<scalar_type, nbredof_> D;
    
    VecE<scalar_type, nbredof_> Fe_ext;

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

template <typename scalar_type, typename real_type, int nbrdof_, int nbredof_>
void ElasticBody<scalar_type, real_type, nbrdof_, nbredof_>::computeForceBalance() {
    const Vec3<scalar_type> aG_local= TG.linear().transpose()*aG;
    const Vec3<scalar_type> omegad_local= TG.linear().transpose()*omegad;
    const Vec3<scalar_type> omega_local= TG.linear().transpose()*omega;

    // M*a
    R= Fext;
    R-= mass * aG_local;
    R-= omegad_local.cross(md);
    for(size_t iedof= 0; iedof < nbredof; iedof++)
        R-= Ct[iedof]*eqdd[iedof];

    MG= Mext;
    MG-= PhiG*omegad_local;
    MG-= md.cross(aG_local);
    for(size_t iedof= 0; iedof < nbredof; iedof++)
        MG-= Cr[iedof]*eqdd[iedof];

    for(size_t iedof= 0; iedof<nbredof; iedof++) {
    	Re[iedof]= Fe_ext[iedof];
        scalar_type prod= Ct[iedof].transpose()*aG_local;
        Re[iedof]-= prod;
        prod= Cr[iedof].transpose()*omegad_local;
        Re[iedof]-= prod;
        prod= Me.row(iedof)*eqdd;
        Re[iedof]-= prod;
    }

    // h volume forces
    R-= omega_local.cross(md.cross(omega_local));
    Vec3<scalar_type> Ct_= Vec3<scalar_type>::Zero();
    for(size_t iedof= 0; iedof<nbredof; iedof++)
        Ct_+= Ct[iedof]*eqd[iedof];
    R-= 2*omega_local.cross(Ct_);

    MG-= omega_local.cross(PhiG*omega_local);
    Mat3<scalar_type> Gr_= Mat3<scalar_type>::Zero();
    for(size_t iedof= 0; iedof<nbredof; iedof++)
        Gr_+= Gr[iedof]*eqd[iedof];
    MG-= Gr_*omega_local;

    Vec6<scalar_type> w;
    w[0]= omega_local[0]*omega_local[0];
    w[1]= omega_local[1]*omega_local[1];
    w[2]= omega_local[2]*omega_local[2];
    w[3]= omega_local[0]*omega_local[1];
    w[4]= omega_local[1]*omega_local[2];
    w[5]= omega_local[0]*omega_local[2];

    Re-= Oe * w;
    for(size_t iedof= 0; iedof<nbredof; iedof++) {
        Vec3<scalar_type> Ge_= Vec3<scalar_type>::Zero();
        for(size_t jedof= 0; jedof<nbredof; jedof++)
            Ge_+= Ge[iedof+jedof*nbredof]*eqd[jedof];
        scalar_type prod= Ge_.transpose()*omega_local;
        Re[iedof]-= prod;
    }

    // h inner forces
	Re.array()-= K.array() * eq.array();
	Re.array()-= D.array() * eqd.array();

    // transform into inertial system
    R= TG.linear() * R;
    MG= TG.linear() * MG;
}

template <typename scalar_type, typename real_type, int nbrdof_, int nbredof_>
void ElasticBody<scalar_type, real_type, nbrdof_, nbredof_>::generalizedBalance(VecQ<scalar_type, nbrdof_> &f) {
	Body<scalar_type, real_type, nbrdof_>::generalizedBalance(f);

	for(size_t iedof= 0; iedof < nbredof; iedof++)
		f[edof[iedof]]-= Re[iedof];
}

template <typename scalar_type, typename real_type, int nbrdof_, int nbredof_>
void ElasticBody<scalar_type, real_type, nbrdof_, nbredof_>::computeGravity(const Vec3<real_type> &g) {
	R+= mass * g;

	Vec3<scalar_type> grav_local= TG.linear().transpose() * g;
	MG+= TG.linear() * md.cross(grav_local);
	for(size_t iedof=0; iedof < nbredof; iedof++) {
        scalar_type prod= Ct[iedof].transpose() * grav_local;
		Re[iedof]+= prod;
    }
}

template <typename scalar_type, typename real_type, int nbrdof_, int nbredof_>
void ElasticBody<scalar_type, real_type, nbrdof_, nbredof_>::setEDOF(const VecQ<scalar_type, nbrdof_> &q, const VecQ<scalar_type, nbrdof_> &qd, const VecQ<scalar_type, nbrdof_> &qdd) {
	for(size_t i= 0; i<nbredof; i++) {
		eq(i)= q(edof[i]);
		eqd(i)= qd(edof[i]);
		eqdd(i)= qdd(edof[i]);
	}
}

#endif /* MULTIBODY_HPP_ */
