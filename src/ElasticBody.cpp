/*
 * ElasticBody.cpp
 *
 *  Created on: 09.03.2017
 *      Author: jgeisler
 */

#include "ElasticBody.hpp"

void ElasticBody::computeForceBalance() {
    const Vec aG_local= TG.linear().transpose()*aG;
    const Vec omegad_local= TG.linear().transpose()*omegad;
    const Vec omega_local= TG.linear().transpose()*omega;

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
    	Re[iedof]= Ree[iedof];
        Re[iedof]-= Ct[iedof].transpose()*aG_local;
        Re[iedof]-= Cr[iedof].transpose()*omegad_local;
        Re[iedof]-= Me.row(iedof)*eqdd;
    }

    // h volume forces
    R-= omega_local.cross(md.cross(omega_local));
    Vec Ct_= Vec::Zero();
    for(size_t iedof= 0; iedof<nbredof; iedof++)
        Ct_+= Ct[iedof]*eqd[iedof];
    R-= 2*omega_local.cross(Ct_);

    MG-= omega_local.cross(PhiG*omega_local);
    Mat Gr_= Mat::Zero();
    for(size_t iedof= 0; iedof<nbredof; iedof++)
        Gr_+= Gr[iedof]*eqd[iedof];
    MG-= Gr_*omega_local;

    Vec6 w;
    w[0]= omega_local[0]*omega_local[0];
    w[1]= omega_local[1]*omega_local[1];
    w[2]= omega_local[2]*omega_local[2];
    w[3]= omega_local[0]*omega_local[1];
    w[4]= omega_local[1]*omega_local[2];
    w[5]= omega_local[0]*omega_local[2];

    Re-= Oe * w;
    for(size_t iedof= 0; iedof<nbredof; iedof++) {
        Vec Ge_= Vec::Zero();
        for(size_t jedof= 0; jedof<nbredof; jedof++)
            Ge_+= Ge[iedof+jedof*nbredof]*eqd[jedof];
        Re[iedof]-= Ge_.transpose()*omega_local;
    }

    // h inner forces
	Re.array()-= K.array() * eq.array();
	Re.array()-= D.array() * eqd.array();

    // transform into inertial system
    R= TG.linear() * R;
    MG= TG.linear() * MG;
}

void ElasticBody::generalizedBalance(Eigen::Ref<VecX> f) {
	Body::generalizedBalance(f);

	for(size_t iedof= 0; iedof < nbredof; iedof++)
		f[edof[iedof]]-= Re[iedof];
}

void ElasticBody::computeGravity(const Eigen::Ref<const Vec> &g) {
	R+= mass * g;

	Vec grav_local= TG.linear().transpose() * g;
	MG+= TG.linear() * md.cross(grav_local);
	for(size_t iedof=0; iedof < nbredof; iedof++)
		Re[iedof]+= Ct[iedof].transpose() * grav_local;
}

void ElasticBody::setEDOF(const Eigen::Ref<const VecX> &q, const Eigen::Ref<const VecX> &qd, const Eigen::Ref<const VecX> &qdd) {
	for(size_t i= 0; i<nbredof; i++) {
		eq(i)= q(edof[i]);
		eqd(i)= qd(edof[i]);
		eqdd(i)= qdd(edof[i]);
	}
}

