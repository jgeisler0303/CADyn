/*
 * GNUPlotIntegrationVisitor.cpp
 *
 *  Created on: 08.04.2017
 *      Author: jgeisler
 */

#include "IntegratorGNUPlotVisitor.hpp"


void IntegratorGNUPlotVisitor::start() {
	if(!system) return;

	if(fname.empty()) {
		fname= system->name;
	}

	writeGNUPlotScript();

	out.open((fname + ".res").c_str());
	system->writeStateVariablesHeader(out);
}

void IntegratorGNUPlotVisitor::step() {
	if(!system) return;

	system->writeStateVariables(out);
}

void IntegratorGNUPlotVisitor::finish() {
	out.close();
}

void IntegratorGNUPlotVisitor::writeGNUPlotScript() {
	std::ofstream s((fname + ".plt").c_str());

	s << "reset" << std::endl;
	s << "set xlabel \"Time [s]\"" << std::endl;
	s << "set grid" << std::endl;
	s << "set term postscript eps color \"Times-Roman\" 20" << std::endl;
	s << "set output \"" << fname << "_q.eps\"" << std::endl;
	s << "set ylabel \"displacements\"" << std::endl;
	s << "plot ";
	for(int i= 0; i < system->nbrdof; i++) {
		if(i>0) s << ", ";
		s << "'" << fname << ".res' using 1:" << i*3+2 << " title '" << system->state_name[i] << "' with line";
	}
	s << std::endl;

	s << "set term postscript eps color \"Times-Roman\" 20" << std::endl;
	s << "set output \"" << fname << "_qd.eps\"" << std::endl;
	s << "set ylabel \"velocities\"" << std::endl;
	s << "plot ";
	for(int i= 0; i < system->nbrdof; i++) {
		if(i>0) s << ", ";
		s << "'" << fname << ".res' using 1:" << i*3+3 << " title 'd" << system->state_name[i] << "' with line";
	}
	s << std::endl;

	s << "set term postscript eps color \"Times-Roman\" 20" << std::endl;
	s << "set output \"" << fname << "_qdd.eps\"" << std::endl;
	s << "set ylabel \"accelerations\"" << std::endl;
	s << "plot ";
	for(int i= 0; i < system->nbrdof; i++) {
		if(i>0) s << ", ";
		s << "'" << fname << ".res' using 1:" << i*3+4 << " title 'dd" << system->state_name[i] << "' with line";
	}
	s << std::endl;
}
