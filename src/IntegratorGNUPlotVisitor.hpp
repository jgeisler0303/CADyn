/*
 * GNUPlotIntegrationVisitor.hpp
 *
 *  Created on: 08.04.2017
 *      Author: jgeisler
 */

#ifndef INTEGRATORGNUPLOTVISITOR_HPP_
#define INTEGRATORGNUPLOTVISITOR_HPP_

#include <fstream>
#include <iostream>
#include <string>
#include "ODEOrder2.hpp"


class IntegratorGNUPlotVisitor: public AbstractIntegratorVisitor {
public:
	IntegratorGNUPlotVisitor(): fname(), out() {};
	IntegratorGNUPlotVisitor(std::string fname_): fname(fname_), out() {};

	virtual void start();
	virtual void step();
	virtual void finish();

	void writeGNUPlotScript();

	std::string fname;
	std::ofstream out;
};

#endif /* INTEGRATORGNUPLOTVISITOR_HPP_ */
