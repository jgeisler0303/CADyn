/*
 * MBSystem.hpp
 *
 *  Created on: 07.03.2017
 *      Author: jgeisler
 */

#ifndef MBSYSTEM_HPP_
#define MBSYSTEM_HPP_

#include "MBSystem.hpp"
#include "MBTypes.hpp"
#include "Body.hpp"
#include "ODEOrder2.hpp"

class MBSystem: public ODEOrder2 {
public:
    MBSystem(int nbrbody_, int nbrdof_, int nbrdep_= 0, int nbrin_= 0, const string &aname= "anonymous_mbs", const string &adesc= "No description"):
    	ODEOrder2(nbrdof_, nbrin_, aname, adesc),
        body(),
        p(nbrdep_),
        pd(nbrdep_),
        pdd(nbrdep_)
        {}

    MBSystem(int nbrbody_, int nbrdof_, int nbrin_= 0, const string &aname= "anonymous_mbs", const string &adesc= "No description"):
    	ODEOrder2(nbrdof_, nbrin_, aname, adesc),
        body(),
        p(0),
        pd(0),
        pdd(0)
    	{}

    virtual ~MBSystem();

    virtual void computeMotion()= 0;
    virtual void computeAppliedEfforts()= 0;
    virtual VecX computeResiduals();
    void computeForceBalances();
    void computeGravity();
    void addSpringForce(double K, double L0, int ibodyA, Vec rA, int ibodyB, Vec rB);
    void addDamperForce(double C, int ibodyA, Vec rA, int ibodyB, Vec rB);

    std::vector<Body*> body;
    VecX p, pd, pdd;
    Vec g= Vec(0.0, 0.0, 0.0);

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

#endif /* MBSYSTEM_HPP_ */
