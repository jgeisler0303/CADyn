/*
 * Body.cpp
 *
 *  Created on: 07.03.2017
 *      Author: jgeisler
 */

#include "Body.hpp"
#include <Eigen/Geometry>

void Body::composeMotion() {
    // TODO throw exception
    if(!rel_base) return;

    Vec erel= rel_base->TG.linear() * TG.translation();
    Vec vrel= rel_base->TG.linear() * vG;
    Vec arel= rel_base->TG.linear() * aG;
    Vec wrel= rel_base->TG.linear() * omega;
    Vec wdrel=rel_base->TG.linear() * omegad;

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

void Body::computeForceBalance() {
    R= Fext;
    R-=mass*aG;

    Mat PhiG_local= TG.linear()*PhiG*TG.linear().transpose();
    MG= Mext;
    MG-= PhiG_local * omegad;
    MG-= omega.cross(PhiG_local*omega);
}

void Body::generalizedBalance(Eigen::Ref<VecX> f) {
	f-= vGpartial.transpose() * R;
	f-= omegapartial.transpose() * MG;
}

void Body::computeGravity(const Eigen::Ref<const Vec> &g) {
	R+= mass * g;
}
