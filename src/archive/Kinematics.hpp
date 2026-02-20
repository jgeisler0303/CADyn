/*
 * Kinematics.hpp
 *
 *  Created on: 07.03.2017
 *      Author: jgeisler
 */

#ifndef KINEMATICS_HPP_
#define KINEMATICS_HPP_

#include "MBTypes.hpp"


class Kinematics {
public:
    Kinematics(size_t nbrdof):
        vGpartial(3, nbrdof),
        omegapartial(3, nbrdof)
	{
    	vGpartial.setZero();
    	omegapartial.setZero();
    	Fext.setZero();
    	Mext.setZero();
	}

	virtual ~Kinematics() {}

	Trans TG;
    Vec vG;
    Vec aG;
    Vec omega;
    Vec omegad;
    Vec R;
    Vec MG;
    Vec Fext;
    Vec Mext;

    Eigen::Matrix<double, 3, Eigen::Dynamic> vGpartial;
    Eigen::Matrix<double, 3, Eigen::Dynamic> omegapartial;
};

#endif /* KINEMATICS_HPP_ */
