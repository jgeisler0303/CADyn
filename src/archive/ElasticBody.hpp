/*
 * ElasticBody.h
 *
 *  Created on: 09.03.2017
 *      Author: jgeisler
 */

#ifndef ELASTICBODY_HPP_
#define ELASTICBODY_HPP_

#include "Body.hpp"

typedef std::vector<Vec, Eigen::aligned_allocator<Vec> > vector_Vec;
typedef std::vector<Mat, Eigen::aligned_allocator<Mat> > vector_Mat;
typedef Eigen::Matrix<double, 6, 1> Vec6;
typedef Eigen::Matrix<double, Eigen::Dynamic, 6> Matx6;

class ElasticBody: public Body {
public:
	ElasticBody(size_t nbrdof, size_t nbredof_, const string &name_= "ElasticBody", const string &desc_= "No description"):
		Body(nbrdof, name_, desc_),
		nbredof(nbredof_),
		edof(nbredof_),
		eq(nbredof_),
		eqd(nbredof_),
		eqdd(nbredof_),
		Ct(nbredof_),
		Cr(nbredof_),
		Me(nbredof_, nbredof_),
		Gr(nbredof_),
		Ge(nbredof_*nbredof_),
		Oe(nbredof_, 6),
		Re(nbredof_),
		Ree(nbredof_),
		K(nbredof_),
		D(nbredof_)
	{
		Ree.setZero();
	}
	virtual ~ElasticBody() {}
    virtual void computeGravity(const Eigen::Ref<const Vec> &g);
	virtual void computeForceBalance();
    virtual void generalizedBalance(Eigen::Ref<VecX> f);
    void setEDOF(const Eigen::Ref<const VecX> &q, const Eigen::Ref<const VecX> &qd, const Eigen::Ref<const VecX> &qdd);

	size_t nbredof;
	std::vector<int> edof;
	VecX eq;
	VecX eqd;
    VecX eqdd;
    Vec md;
    vector_Vec Ct;
    vector_Vec Cr;
    MatX Me;

    vector_Mat Gr;
    vector_Vec Ge;
    Matx6 Oe;

    VecX Re;
    VecX Ree;

    VecX K;
    VecX D;

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

#endif /* ELASTICBODY_HPP_ */
