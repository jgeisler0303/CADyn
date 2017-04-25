/*
 * MBSystem.cpp
 *
 *  Created on: 07.03.2017
 *      Author: jgeisler
 */

#include "MBSystem.hpp"

MBSystem::~MBSystem() {
    for(size_t ibody= 0; ibody < body.size(); ibody++)
        delete body[ibody];
}

void MBSystem::computeForceBalances() {
    for (size_t ibody= 0; ibody < body.size(); ibody++)
        body[ibody]->computeForceBalance();
}

void MBSystem::computeGravity() {
    for (size_t ibody= 0; ibody < body.size(); ibody++)
        body[ibody]->computeGravity(g);
}

VecX MBSystem::computeResiduals() {
	VecX f_(nbrdof);

    computeMotion();
    computeAppliedEfforts();
    computeForceBalances();
    computeGravity();

    f_.setZero();
    for (size_t ibody= 0; ibody < body.size(); ibody++)
        body[ibody]->generalizedBalance(f_);

    return f_;
}

void MBSystem::addSpringForce(double K, double L0, int ibodyA, Vec rA, int ibodyB, Vec rB) {
	Vec uAB= body[ibodyB]->TG*rB - body[ibodyA]->TG*rA;
	double L= uAB.norm();
	uAB.normalize();
	Vec force= K*(L-L0)*uAB;
	body[ibodyA]->Fext+= force;
	body[ibodyA]->Mext+= (body[ibodyA]->TG.linear()*rA).cross(force);
	body[ibodyB]->Fext-= force;
	body[ibodyB]->Mext-= (body[ibodyB]->TG.linear()*rB).cross(force);
}

void MBSystem::addDamperForce(double C, int ibodyA, Vec rA, int ibodyB, Vec rB) {
	Vec uAB= body[ibodyB]->TG*rB - body[ibodyA]->TG*rA;
	uAB.norm();
	Vec vrel= body[ibodyB]->vG + body[ibodyB]->omega.cross(body[ibodyB]->TG.linear()*rB) - body[ibodyA]->vG - body[ibodyA]->omega.cross(body[ibodyA]->TG.linear()*rA);
	Vec force= C*(vrel.dot(uAB))*uAB;
	body[ibodyA]->Fext+= force;
	body[ibodyA]->Mext+= (body[ibodyA]->TG.linear()*rA).cross(force);
	body[ibodyB]->Fext-= force;
	body[ibodyB]->Mext-= (body[ibodyB]->TG.linear()*rB).cross(force);
}
