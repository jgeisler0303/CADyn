/*
 * IntegratorCompareVisitor.hpp
 *
 *  Created on: 19.04.2017
 *      Author: jgeisler
 */

#ifndef INTEGRATORCOMPAREVISITOR_HPP_
#define INTEGRATORCOMPAREVISITOR_HPP_

#include "ODEOrder2.hpp"
#include <string>
#include <fstream>
#include <iostream>

class IntegratorCompareVisitor: public AbstractIntegratorVisitor {
public:
	IntegratorCompareVisitor(std::string fname_): fname(fname_) {};

	virtual void start();
	virtual void step();
	virtual void finish();

	double rms_q(int i);
	double rms_qd(int i);
	double rms_qdd(int i);
	double nrms_q(int i);
	double nrms_qd(int i);
	double nrms_qdd(int i);
	double total_nrms_q();
	double total_nrms_qd();
	double total_nrms_qdd();
	double total_nrms();

	void printResult();

	std::string fname;
	std::ifstream ref_file;

	VecX e_q;
	VecX e_qd;
	VecX e_qdd;

	VecX min_q, max_q, min_qd, max_qd, min_qdd, max_qdd;

	int nline;
};

#endif /* INTEGRATORCOMPAREVISITOR_HPP_ */
