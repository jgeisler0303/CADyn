/*
 * Body.h
 *
 *  Created on: 07.03.2017
 *      Author: jgeisler
 */

#ifndef BODY_H_
#define BODY_H_

#include "Kinematics.hpp"


class Body: public Kinematics {
public:
    Body(size_t nbrdof, const string &aname= "Body", const string &adesc= "No description"): Kinematics(nbrdof),
		name(aname),
		description(adesc),
        mass(0.0),
        PhiG()
        {}

    virtual ~Body() { /* delete relative; */ };

    string name;
    string description;

    void composeMotion();
    virtual void computeGravity(const Eigen::Ref<const Vec> &g);
	virtual void computeForceBalance();
    virtual void generalizedBalance(Eigen::Ref<VecX> f);

    void applyForceInLocal(const Eigen::Ref<const Vec> &r, const Eigen::Ref<const Vec> &F);
    void applyForceIn0(const Eigen::Ref<const Vec> &r, const Eigen::Ref<const Vec> &F);
        
    double mass;
    Mat PhiG;

    Kinematics *rel_base= nullptr;
    // Kinematics *relative= nullptr;

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

#endif /* BODY_H_ */
