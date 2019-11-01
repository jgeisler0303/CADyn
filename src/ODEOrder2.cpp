/*
 * ODEOrder2.cpp
 *
 *  Created on: 10.03.2017
 *      Author: jgeisler
 */

#include "ODEOrder2.hpp"
#include <cmath>
#include <algorithm>
#include <sstream>

ODEOrder2::ODEOrder2(int nbrdof_, int nbrin_, const string &aname, const string &adesc):
        name(aname),
        description(adesc),
		nbrdof(nbrdof_),
		nbrin(nbrin_),
		state_name(nbrdof_),
		in_name(nbrin_),
		Jacobian(nbrdof_, nbrdof_),
		LU(),
        q(nbrdof_),
        qd(nbrdof_),
        qdd(nbrdof_),
		f(nbrdof_),
		doflocked(nbrdof_, false),
		u(nbrin_),
		t(0.0)
	{
		q.setZero();
		qd.setZero();
		qdd.setZero();
		u.setZero();
		for(int i= 0; i<nbrdof; i++) {
			std::stringstream ss;
			ss << "q_" << i;
			state_name[i]= ss.str();
		}
	}

void ODEOrder2::writeStateVariablesHeader(std::ostream &OutFile) {
	OutFile << " time ";

	for(int idof= 0; idof < nbrdof; idof++)
		OutFile << " q" << idof << " qd" << idof << " qdd" << idof;

	for(int iinput= 0; iinput < nbrin; iinput++)
		OutFile << " u" << iinput;

	OutFile << std::endl;
}

void ODEOrder2::writeStateVariables(std::ostream &OutFile) {
	OutFile << " " << t;

	for(int idof= 0;idof < nbrdof; idof++)
		OutFile << " " << q[idof] << " " << qd[idof] << " " << qdd[idof];

	for(int iinput= 0; iinput < nbrin; iinput++)
		OutFile << " " << u[iinput];

	OutFile << std::endl;
}

VecX ODEOrder2::computeResidualsInt() {
	VecX f= computeResiduals();
	for(int idof= 0; idof < nbrdof; idof++)
		if(doflocked[idof]) f[idof]= 0.0;

	return f;
}

void ODEOrder2::calcJacobian(double alphaM, double alphaC, double alphaK, double tol) {
	VecX foff;
	f= computeResidualsInt();

	for(int jddl= 0; jddl < nbrdof; jddl++) {
		double q_= q[jddl];
		double qd_= qd[jddl];
		double qdd_= qdd[jddl];

		q[jddl]+= alphaK*tol;
		qd[jddl]+= alphaC*tol;
		qdd[jddl]+= alphaM*tol;

		if(!doflocked[jddl]) {
			foff= computeResidualsInt();
			Jacobian.col(jddl)= (foff - f)/tol;
		} else
			Jacobian.col(jddl).setZero();

		q[jddl]= q_;
		qd[jddl]= qd_;
		qdd[jddl]= qdd_;
	}
}

bool ODEOrder2::staticEquilibrium() {
	double err;

	int nstep= 0;
	do {
		nstep++;
		if((nstep%10)==1) {
			calcJacobian(100.0, 0.0, 1.0, 1.0E-5);
			LU.compute(Jacobian);
		} else
			f= computeResidualsInt();

		VecX qdd_corr= LU.solve(f);
		q-= qdd_corr;
		err= qdd_corr.norm();
	}
	while((nstep<1000) && (err > (1E-8 * sqrt(1.0*nbrdof))));

	return(nstep<=1000);
}

int ODEOrder2::newmarkOneStep(double h, double &errq, bool hmodified) {
	const double Beta= 0.25;
	const double Gamma=0.5;

	VecX qdd_sto= qdd;
	t+= h;
	q+= h*qd + 0.5*h*h*qdd;
	qd+= h*qdd;

	int istepjac= 0;
	if(hmodified) istepjac= 1;

	int nstep=0;
	double err;
	do {
		nstep++;
		if((nstep%4)==istepjac) {
			calcJacobian(1.0, Gamma*h, Beta*h*h, 1E-2);
			LU.compute(Jacobian);
		} else
			f= computeResidualsInt();

		VecX qdd_corr= LU.solve(f);
		err= qdd_corr.norm() / (sqrt(1.0*nbrdof) * (1.0 + qdd.norm()));

		qdd-= qdd_corr;
		qdd_corr*= Gamma*h;
		qd-= qdd_corr;
		qdd_corr*= Beta*h/Gamma;
		q-= qdd_corr;
	} while((nstep<10) && (err>1E-8));

	if(nstep==10) return 1;

	qdd_sto-= qdd;
	errq= 0.0;
	for(int iddl=0; iddl<nbrdof; iddl++) {
		double errtmp, absqiddl;
		errtmp= h*h*fabs(qdd_sto[iddl]) / 12.0;
		absqiddl= fabs(q[iddl]);
		absqiddl= std::max(absqiddl, h*fabs(qd[iddl]));
		absqiddl= std::max(absqiddl, h*h*fabs(qdd[iddl]));
		if((absqiddl*RelTol)>AbsTol) errtmp= errtmp*AbsTol/(absqiddl*RelTol);
		errq+= errtmp*errtmp;
	}
	errq= sqrt(errq/nbrdof);

	if (!std::isfinite(errq)) return 2;
	return 0;
}

bool ODEOrder2::newmarkInterval(double tfinal, double &h, double hmax) {
	bool hchanged= true;
	double errq;
    n_steps= 0;
    n_back_steps= 0;

	if (h>hmax) h= hmax;
	while(t < tfinal) {
		if((t+1.4*h) >= tfinal) {
			h= tfinal-t;
			hchanged= true;
		}
		double timesto= t;

		VecX q_sto= q;
		VecX qd_sto= qd;
		VecX qdd_sto= qdd;

		int code= newmarkOneStep(h, errq, hchanged);
        n_steps++;
		hchanged= false;

		if((code) || (errq>AbsTol))	{
			if((code==1) || (code==2) || !std::isfinite(errq))
				h*= 0.25;
			else
				h*= sqrt((0.21*AbsTol + 0.04*errq) / errq);
			hchanged= true;
			t= timesto;

			q= q_sto;
			qd= qd_sto;
			qdd= qdd_sto;
			if(h<hminmin)
				return false;
            
            n_back_steps++;
		} else {
			if((errq < (0.1*AbsTol)) && (h<hmax)) {
				h*= sqrt(AbsTol/(2.1*errq + 0.04*AbsTol));
				if(h>hmax) h= hmax;
				hchanged= true;
			}
		}
	}
	return true;
}

bool ODEOrder2::newmarkIntegration(double tfinal, double hsave, double hmax, AbstractIntegratorVisitor *visitor) {
	t= 0.0;
	double h= 0.0;
	double dummy;

	if(visitor) {
		visitor->setSystem(this);
		visitor->start();
	}

	newmarkOneStep(0.0, dummy);
	if(visitor) visitor->step();

	int ipas= 0;
	bool res= true;
	h= hmax;
	while(t<tfinal)	{
		ipas++;
		if(!newmarkInterval(hsave*ipas, h, hmax)) {
			res= false;
			break;
		}
		if(visitor) visitor->step();
	}

	if(visitor) visitor->finish();
	return res;
}


//void SaveLinearizedSystem(int *doflocked)
//
//// Saves mass, damping and stiffness matrices and eventually
//// inputs contribution
//
//{
//  int idof,jdof;
//  char *nom;
//  nom=new char[strlen(application)+5];
//
//  cout << "Constructing and saving stiffness matrix ...";
//  Calcule_Matrice_J(0,0,1,1E-5,doflocked);
//  strcpy(nom,application);
//  strcat(nom,".kk");
//  ofstream stifffile(nom);
//  if (!stifffile) EasyDynsimError("Cannot create stiffness matrix file");
//  for (idof=0; idof<nbrdof; idof++)
//     {
//     for (jdof=0;jdof<nbrdof; jdof++)
//	 stifffile << " " << gsl_matrix_get(jac,idof,jdof);
//     stifffile << endl;
//     }
//  stifffile.close();
//  cout << " Done\n";
//
//  cout << "Constructing and saving damping matrix ...";
//  Calcule_Matrice_J(0,1,0,1E-3,doflocked);
//  strcpy(nom,application);
//  strcat(nom,".cc");
//  ofstream dampfile(nom);
//  if (!dampfile) EasyDynsimError("Cannot create damping matrix file");
//  for (idof=0; idof<nbrdof; idof++)
//     {
//     for (jdof=0;jdof<nbrdof; jdof++)
//	 dampfile << " " << gsl_matrix_get(jac,idof,jdof);
//     dampfile << endl;
//     }
//  dampfile.close();
//  cout << " Done\n";
//
//  cout << "Constructing and saving mass matrix ...";
//  Calcule_Matrice_J(1,0,0,1,doflocked);
//  strcpy(nom,application);
//  strcat(nom,".mm");
//  ofstream massfile(nom);
//  if (!massfile) EasyDynsimError("Cannot create mass matrix file");
//  for (idof=0; idof<nbrdof; idof++)
//     {
//     for (jdof=0;jdof<nbrdof; jdof++)
//	 massfile << " " << gsl_matrix_get(jac,idof,jdof);
//     massfile << endl;
//     }
//  massfile.close();
//  cout << " Done\n";
//
//  if (nbrinput>0)
//    {
//    cout << "Constructing and saving input contribution matrix ...";
//    double **FF, tol=1e-3;
//    int iinput;
//    FF=new double*[nbrdof];
//    for (idof=0; idof<nbrdof;idof++) FF[idof]=new double[nbrinput];
//    for (iinput=0; iinput<nbrinput; iinput++)
//      {
//      u[iinput]+=tol;
//      ComputeResidual();
//      for (idof=0; idof<nbrdof;idof++) FF[idof][iinput]=0.5*f[idof]/tol;
//      u[iinput]-=2*tol;
//      ComputeResidual();
//      for (idof=0; idof<nbrdof;idof++) FF[idof][iinput]-=0.5*f[idof]/tol;
//      u[iinput]+=tol;
//      }
//    strcpy(nom,application);
//    strcat(nom,".ff");
//    ofstream fffile(nom);
//    if (!fffile) EasyDynsimError("Cannot create influence matrix file");
//    for (idof=0;idof<nbrdof; idof++)
//      {
//      for (iinput=0; iinput<nbrinput; iinput++)
//        fffile << " " << -FF[idof][iinput];
//      fffile << endl;
//      }
//    fffile.close();
//    cout << " Done\n";
//    for (idof=0;idof<nbrdof; idof++) delete FF[idof];
//    delete FF;
//    }
//}

//void ComputePoles(int *doflocked, double freqmin, double freqmax)
//// Procede au calcul des poles pour la configuration actuelle
//
//{
//  int NN,INFO,*Liste;
//  gsl_matrix *AA,*BB;
//  gsl_matrix_complex *VEC;
//  gsl_vector_complex *ALPHA;
//  gsl_vector *BETA;
//  gsl_eigen_genv_workspace *WORK;
//  NN=2*nbrdof;
//  // Memory allocation
//  AA=gsl_matrix_alloc(NN,NN);
//  BB=gsl_matrix_alloc(NN,NN);
//  VEC=gsl_matrix_complex_alloc(NN,NN);
//  ALPHA=gsl_vector_complex_alloc(NN);
//  BETA=gsl_vector_alloc(NN);
//  WORK=gsl_eigen_genv_alloc(NN);
//  Liste=new int[NN];
//
//  // Building eigenvalue problem matrices
//  int iddl,jddl,ipole,jpole;
//  // Zeroing matrices
//  gsl_matrix_set_zero(AA);
//  gsl_matrix_set_zero(BB);
//  // Adding stiffness matrix contribution
//  cout << "Constructing Stiffness matrix ...";
//  Calcule_Matrice_J(0,0,1,1E-5,doflocked);
//  cout << " Done\n";
//  for (iddl=0; iddl<nbrdof; iddl++) for (jddl=0;jddl<nbrdof; jddl++)
//     gsl_matrix_set(AA,iddl,jddl,-gsl_matrix_get(jac,iddl,jddl));
//
//  // Adding damping matrix contribution
//  cout << "Constructing Damping matrix ...";
//  Calcule_Matrice_J(0,1,0,1E-3,doflocked);
//  cout << " Done\n";
//  for (iddl=0; iddl<nbrdof; iddl++) for (jddl=0;jddl<nbrdof; jddl++)
//     gsl_matrix_set(BB,iddl,jddl,gsl_matrix_get(jac,iddl,jddl));
//
//  // Adding mass matrix contribution
//  cout << "Constructing Mass matrix ...";
//  Calcule_Matrice_J(1,0,0,1,doflocked);
//  cout << " Done\n";
//  for (iddl=0; iddl<nbrdof; iddl++) for (jddl=0;jddl<nbrdof; jddl++)
//     {
//     gsl_matrix_set(AA,nbrdof+iddl,nbrdof+jddl,gsl_matrix_get(jac,iddl,jddl));
//     gsl_matrix_set(BB,nbrdof+iddl,jddl,gsl_matrix_get(jac,iddl,jddl));
//     gsl_matrix_set(BB,iddl,nbrdof+jddl,gsl_matrix_get(jac,iddl,jddl));
//     }
//
//  // Calling GSL routine for nonsymmetric problem
//  INFO=gsl_eigen_genv(AA,BB,ALPHA,BETA,VEC,WORK);
//  gsl_eigen_genv_free(WORK); // immediately freeing work memory
//  if (INFO!=0) // did everything happen properly ?
//   {
//   cout << "Code returned by gsl_eigen_genv: " << INFO << endl;
//   EasyDynsimError("gsl_eigen_genv failed");
//   }
//
//  // Reporting eigen values to ALPHA
//  double beta; gsl_complex alpha;
//  for (ipole=0; ipole<NN; ipole++)
//    {
//    beta=gsl_vector_get(BETA,ipole);
//    alpha=gsl_vector_complex_get(ALPHA,ipole);
//    if (fabs(beta)<1E-15)
//      {
//      beta=1.0;
//      alpha=gsl_complex_rect(1E30,1E30);
//      }
//    alpha=gsl_complex_div_real(alpha,beta);
//    gsl_vector_set(BETA,ipole,beta);
//    gsl_vector_complex_set(ALPHA,ipole,alpha);
//    }
//  // Building ordered list of poles
//  double sigmai,sigmaj, omegai, omegaj;
//  for (ipole=0; ipole<NN; ipole++) Liste[ipole]=ipole;
//  for (ipole=0; ipole<NN-1; ipole++) for (jpole=ipole; jpole<NN; jpole++)
//    {
//    alpha=gsl_vector_complex_get(ALPHA,Liste[ipole]);
//    omegai=GSL_IMAG(alpha);
//    sigmai=GSL_REAL(alpha);
//    alpha=gsl_vector_complex_get(ALPHA,Liste[jpole]);
//    omegaj=GSL_IMAG(alpha);
//    sigmaj=GSL_REAL(alpha);
//    if (omegaj<omegai)
//           {
//           int itmp=Liste[ipole];
//           Liste[ipole]=Liste[jpole];
//           Liste[jpole]=itmp;
//           }
//    if (fabs(omegaj-omegai)<1E-15)
//      if (sigmaj<sigmai)
//           {
//           int itmp=Liste[ipole];
//           Liste[ipole]=Liste[jpole];
//           Liste[jpole]=itmp;
//           }
//    }
//  // printing eigenvalues if in debug mode
//  if (DEBUG) for (ipole=0; ipole<NN; ipole++)
//    {
//    alpha=gsl_vector_complex_get(ALPHA,Liste[ipole]);
//    omegai=GSL_IMAG(alpha);
//    sigmai=GSL_REAL(alpha);
//    if (omegai>-1E-15)
//        cout << "Pole " << Liste[ipole] << "=" << sigmai
//             << "+-i" << omegai << "\n";
//    }
//  // Building complete results file
//  int numpole, FREQMINMAX=0;
//  double polenorm,freq,amor,coeff;
//  char *nom;
//  if ((freqmin!=0) || (freqmax!=1E30)) FREQMINMAX=1;
//
//  nom=new char[strlen(application)+5];
//  strcpy(nom,application);
//  strcat(nom,".lst"); // file with eigenvalues
//  ofstream lstfile(nom);
//  if (!lstfile)
//    {
//    cout<<"\n Fichier "<<nom<<" impossible a ouvrir";
//    exit(1);
//    }
//  lstfile << "  POLE       ALPHA         OMEGA        "
//        << "FREQ(Hz)   DAMP_RATIO \n";
//
//  ofstream ls2file;
//  if (FREQMINMAX)
//    {
//    strcpy(nom,application);
//    strcat(nom,".ls2");
//    ls2file.open(nom);
//    if (!ls2file)
//      {
//      cout<<"\n Cannot open " << nom;
//      exit(1);
//      }
//    }
//  strcpy(nom,application);
//  strcat(nom,".mod");
//  ofstream modfile(nom);
//  if (!modfile)
//    {
//    cout<<"\n Cannot open " << nom;
//    exit(1);
//    }
//
//  for(ipole=0;ipole<NN;ipole++)
//    {
//    numpole=Liste[ipole];
//    alpha=gsl_vector_complex_get(ALPHA,numpole);
//    omegai=GSL_IMAG(alpha);
//    sigmai=GSL_REAL(alpha);
//    polenorm=hypot(sigmai,omegai);
//    freq=omegai/(2*3.14159265);
//    if (polenorm<1E-15) polenorm=1;
//    amor=-sigmai/polenorm;
//    if (omegai>-1E-15) lstfile << setw(5) << ipole+1
//            << " " << setw(13) << sigmai
//            << " " << setw(13) << omegai
//            << " " << setw(13) << freq
//            << " " << setw(13) << amor <<"\n";
//    if ((freq>=freqmin) && (freq<=freqmax))
//      {
//      if (FREQMINMAX) ls2file << setw(5) << ipole+1
//              << " " << setw(13) << sigmai
//              << " " << setw(13) << omegai
//              << " " << setw(13) << freq
//              << " " << setw(13) << amor <<"\n";
//      modfile << setw(5) << ipole+1
//              << " " << setw(13) << sigmai
//              << " " << setw(13) << omegai
//              << " " << setw(13) << freq
//              << " " << setw(13) << amor <<"\n";
//      // Saving mode
//      if (omegai>-1E-15) // conjugate poles are saved only once
//        {
//        for (iddl=0; iddl<nbrdof; iddl++)
//	  {
//	  alpha=gsl_matrix_complex_get(VEC,nbrdof+iddl,numpole);
//          modfile << GSL_REAL(alpha) << " " << GSL_IMAG(alpha) << endl;
//          }
//        }
//      }
//    }
//  lstfile.close();
//  if (FREQMINMAX) ls2file.close();
//  modfile.close();
//
//  // Cleaning memeory
//  gsl_matrix_free(AA);
//  gsl_matrix_free(BB);
//  gsl_matrix_complex_free(VEC);
//  gsl_vector_complex_free(ALPHA);
//  gsl_vector_free(BETA);
//  delete Liste;
//  }
