/*
 * IntegratorCompareVisitor.cpp
 *
 *  Created on: 19.04.2017
 *      Author: jgeisler
 */

#include "IntegratorCompareVisitor.hpp"
#include <sstream>
#include <stdexcept>

void IntegratorCompareVisitor::start() {
	if(!system) return;

	ref_file.open(fname.c_str());


	std::string line;
	std::getline(ref_file, line);
	if(ref_file.fail()) throw std::domain_error("Reference file is empty.");

	std::istringstream in(line);

	std::string time_str;
	in >> time_str;
	if(in.fail() || time_str.compare("time")!=0) {
		std::stringstream err_msg;
		err_msg << "Invalid header in reference file. 'time' not found. Found: " << time_str;
		throw std::domain_error(err_msg.str());
	}
	for(int idof= 0; idof < system->nbrdof; idof++) {
		std::string q_name;
		std::string qd_name;
		std::string qdd_name;

		in >> q_name >> qd_name >> qdd_name;

		if(in.fail()) {
			std::stringstream err_msg;
			err_msg << "Invalid header in reference file. Names of state " << idof <<" not found.";
			throw std::domain_error(err_msg.str());
		}
		if(q_name.compare(0, 1, "q")!=0) {
			std::stringstream err_msg;
			err_msg << "Invalid header in reference file. Names of state " << idof <<" don't start with 'q'. Found: " << q_name;
			throw std::domain_error(err_msg.str());
		}
		if(qd_name.compare(0, 2, "qd")!=0) {
			std::stringstream err_msg;
			err_msg << "Invalid header in reference file. Names of state " << idof <<" don't start with 'qd'. Found: " << qd_name;
			throw std::domain_error(err_msg.str());
		}
		if(qdd_name.compare(0, 3, "qdd")!=0) {
			std::stringstream err_msg;
			err_msg << "Invalid header in reference file. Names of state " << idof <<" don't start with 'qdd'. Found: " << qdd_name;
			throw std::domain_error(err_msg.str());
		}
	}

	e_q.setZero(system->nbrdof);
	e_qd.setZero(system->nbrdof);
	e_qdd.setZero(system->nbrdof);

	min_q.setZero(system->nbrdof);
	max_q.setZero(system->nbrdof);
	min_qd.setZero(system->nbrdof);
	max_qd.setZero(system->nbrdof);
	min_qdd.setZero(system->nbrdof);
	max_qdd.setZero(system->nbrdof);

	nline= 0;
}

void IntegratorCompareVisitor::step() {
	nline++;

	std::string line;
	std::getline(ref_file, line);
	if(ref_file.fail()){
		std::stringstream err_msg;
		err_msg << "Reference file too short at line " << nline;
		throw std::domain_error(err_msg.str());
	}

	std::istringstream in(line);
	double time;
	in >> time;
	if(in.fail()) {
		std::stringstream err_msg;
		err_msg << "Invalid data in reference file. Values for time not found at line " << nline << ".";
		throw std::domain_error(err_msg.str());
	}
	if(fabs(time-system->t)>1e-6) {
		std::stringstream err_msg;
		err_msg << "Invalid data in reference file. Values for time diverge at line " << nline << ". Should be " << system->t << " is " << time;
		throw std::domain_error(err_msg.str());
	}

	for(int idof= 0; idof < system->nbrdof; idof++) {
		double ref_q, ref_qd, ref_qdd;

		in >> ref_q >> ref_qd >> ref_qdd;

		if(in.fail()) {
			std::stringstream err_msg;
			err_msg << "Invalid data in reference file. Values for state " << idof <<" not found at line " << nline << ".";
			throw std::domain_error(err_msg.str());
		}

		e_q(idof)+= pow(ref_q-system->q(idof), 2);
		e_qd(idof)+= pow(ref_qd-system->qd(idof), 2);
		e_qdd(idof)+= pow(ref_qdd-system->qdd(idof), 2);

		min_q(idof)= std::min(min_q(idof), ref_q);
		max_q(idof)= std::max(max_q(idof), ref_q);
		min_qd(idof)= std::min(min_qd(idof), ref_qd);
		max_qd(idof)= std::max(max_qd(idof), ref_qd);
		min_qdd(idof)= std::min(min_qdd(idof), ref_qdd);
		max_qdd(idof)= std::max(max_qdd(idof), ref_qdd);
	}
}

void IntegratorCompareVisitor::finish() {
}

double IntegratorCompareVisitor::rms_q(int i) {
	if(nline==0)
		return 0.0;
	else
		return sqrt(e_q(i))/double(nline);
}

double IntegratorCompareVisitor::rms_qd(int i) {
	if(nline==0)
		return 0.0;
	else
		return sqrt(e_qd(i))/double(nline);
}

double IntegratorCompareVisitor::rms_qdd(int i) {
	if(nline==0)
		return 0.0;
	else
		return sqrt(e_qdd(i))/double(nline);
}

double IntegratorCompareVisitor::nrms_q(int i) {
	double range= max_q(i)-min_q(i);
	if(range<1e-8)
		return rms_q(i);
	else
		return rms_q(i)/range;
}

double IntegratorCompareVisitor::nrms_qd(int i) {
	double range= max_qd(i)-min_qd(i);
	if(range<1e-8)
		return rms_qd(i);
	else
		return rms_qd(i)/range;
}

double IntegratorCompareVisitor::nrms_qdd(int i) {
	double range= max_qdd(i)-min_qdd(i);
	if(range<1e-8)
		return rms_qdd(i);
	else
		return rms_qdd(i)/range;
}

double IntegratorCompareVisitor::total_nrms_q() {
	if(e_q.size()<=0) return 0.0;

	double sum= 0.0;
	for(int i= 0; i<e_q.size(); i++)
		sum+= nrms_q(i);

	return sum/double(e_q.size());
}

double IntegratorCompareVisitor::total_nrms_qd() {
	if(e_qd.size()<=0) return 0.0;

	double sum= 0.0;
	for(int i= 0; i<e_qd.size(); i++)
		sum+= nrms_qd(i);

	return sum/double(e_qd.size());
}

double IntegratorCompareVisitor::total_nrms_qdd() {
	if(e_qdd.size()<=0) return 0.0;

	double sum= 0.0;
	for(int i= 0; i<e_qdd.size(); i++)
		sum+= nrms_qdd(i);

	return sum/double(e_qdd.size());
}

double IntegratorCompareVisitor::total_nrms() {
	return (total_nrms_q() + total_nrms_qd() + total_nrms_qdd()) / 3.0;
}

void IntegratorCompareVisitor::printResult() {
	std::cout << "Absolute errors of q: ";
	for(int i= 0; i<e_q.size(); i++)
		std::cout << rms_q(i) << ", ";
	std::cout << std::endl;

	std::cout << "Normalized errors of q: ";
	for(int i= 0; i<e_q.size(); i++)
        std::cout << nrms_q(i) << ", ";
	std::cout << std::endl;

	std::cout << "Absolute errors of qd: ";
	for(int i= 0; i<e_qd.size(); i++)
        std::cout << rms_qd(i) << ", ";
	std::cout << std::endl;

	std::cout << "Normalized errors of qd: ";
	for(int i= 0; i<e_qd.size(); i++)
        std::cout << nrms_qd(i) << ", ";
	std::cout << std::endl;

	std::cout << "Absolute errors of qdd: ";
	for(int i= 0; i<e_qd.size(); i++)
        std::cout << rms_qdd(i) << ", ";
	std::cout << std::endl;

	std::cout << "Normalized errors of qdd: ";
	for(int i= 0; i<e_qdd.size(); i++)
        std::cout << nrms_qdd(i) << ", ";
	std::cout << std::endl;

	std::cout << "Total normalized errors. q: " << total_nrms_q() << ", qd: " << total_nrms_qd() << ", qdd: " << total_nrms_qdd() << std::endl;
	std::cout << "Total normalized error: " << total_nrms() << std::endl << std::endl;
}
