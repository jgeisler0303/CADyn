#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define makeMBSystemClass(c) c ## System
#define exp_makeMBSystemClass(c) makeMBSystemClass(c)
#define MBSystemClass exp_makeMBSystemClass(MBSystem)

#define stringify(x)  #x
#define expand_and_stringify(x) stringify(x ## System2.hpp)
#define INCL_FILE_STR(x) expand_and_stringify(x)
#include INCL_FILE_STR(MBSystem)

#include "mex.h"
#ifndef  HAVE_OCTAVE
#include "matrix.h"
#endif

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    if(nrhs!=5) { mexErrMsgIdAndTxt("CADyn:InvalidArgument", "Wrong number of arguments. Expecting (x0, dx0, ddx0, u, p)"); return; }
    if(nlhs<1 || nlhs>4) { mexErrMsgIdAndTxt("CADyn:InvalidArgument", "Wrong number of return values. Expecting [f, {M, {C, {K}}}}]"); return; }
    
    if(!mxIsDouble(prhs[0]) || mxGetNumberOfElements(prhs[0])!=MBSystemClass::nStates) { mexErrMsgIdAndTxt("CADyn:InvalidArgument", "Wrong number of elements in 'x0' (%d expected)", MBSystemClass::nStates); return; }
    if(!mxIsDouble(prhs[1]) || mxGetNumberOfElements(prhs[1])!=MBSystemClass::nStates) { mexErrMsgIdAndTxt("CADyn:InvalidArgument", "Wrong number of elements in 'dx0' (%d expected)", MBSystemClass::nStates); return; }
    if(!mxIsDouble(prhs[2]) || mxGetNumberOfElements(prhs[2])!=MBSystemClass::nStates) { mexErrMsgIdAndTxt("CADyn:InvalidArgument", "Wrong number of elements in 'ddx0' (%d expected)", MBSystemClass::nStates); return; }
    if(!mxIsDouble(prhs[3]) || mxGetNumberOfElements(prhs[3])!=MBSystemClass::nInputs) { mexErrMsgIdAndTxt("CADyn:InvalidArgument", "Wrong number of elements in 'u' (%d expected)", MBSystemClass::nInputs); return; }
    
    const mxArray *mxParams= prhs[4];
    if(!mxIsStruct(mxParams) || mxGetNumberOfElements(mxParams)!=1) {
        mexErrMsgIdAndTxt("CADyn:InvalidArgument", "Input p must be a scalar struct.\n");
        return;
    }
    
    MBSystemClass system;
    
    for(const auto &iter : system.param.info_map) {
        const mxArray *mxParam;
        if((mxParam= mxGetField(mxParams, 0, iter.first.c_str()))==NULL) {
            mexErrMsgIdAndTxt("CADyn:InvalidArgument", "Required parameter '%s' is not member of parameters struct.\n", iter.first.c_str());
            return;
        }
        int m_= mxGetM(mxParam);
        int n_= mxGetN(mxParam);
        if(mxIsSparse(mxParam) || !mxIsDouble(mxParam) || m_!=iter.second.nrows || n_!=iter.second.ncols) {
            mexErrMsgIdAndTxt("CADyn:InvalidArgument", "Parameter name '%s' must be a %dx%d matrix.\n", iter.first.c_str(), iter.second.nrows, iter.second.ncols);
            return;
        }
        system.param.setParam(iter.first.c_str(), mxGetPr(mxParam));
    }
    system.precalcConsts();
    
    double *x0= mxGetPr(prhs[0]);
    for(int i=0; i<MBSystemClass::nStates; ++i)
        system.q(i)= x0[i];
    
    double *dx0= mxGetPr(prhs[1]);
    for(int i=0; i<MBSystemClass::nStates; ++i)
        system.qd(i)= dx0[i];
    
    double *ddx0= mxGetPr(prhs[2]);
    for(int i=0; i<MBSystemClass::nStates; ++i)
        system.qdd(i)= ddx0[i];
    
    double *u_sim= mxGetPr(prhs[3]);;
    for(int i=0; i<MBSystemClass::nInputs; ++i)
        system.u(i)= u_sim[i];
    
    system.calcJacobian(1.0, 1.0, 1.0);
    
    
    plhs[0]= mxCreateDoubleMatrix(MBSystemClass::nStates, 1, mxREAL);
    for(int i=0; i<MBSystemClass::nStates; ++i)
        mxGetPr(plhs[0])[i]= system.f(i);
        
    if(nlhs>1) {
        plhs[1]= mxCreateDoubleMatrix(MBSystemClass::nStates, MBSystemClass::nStates, mxREAL);
        for(int irow=0; irow<MBSystemClass::nStates; ++irow)
            for(int icol=0; icol<MBSystemClass::nStates; ++icol)
                mxGetPr(plhs[1])[irow + icol*MBSystemClass::nStates]= system.M(irow, icol);
    }
    
    if(nlhs>2) {
        plhs[2]= mxCreateDoubleMatrix(MBSystemClass::nStates, MBSystemClass::nStates, mxREAL);
        for(int irow=0; irow<MBSystemClass::nStates; ++irow)
            for(int icol=0; icol<MBSystemClass::nStates; ++icol)
                mxGetPr(plhs[2])[irow + icol*MBSystemClass::nStates]= system.C(irow, icol);
    }
    
    if(nlhs>3) {
        plhs[3]= mxCreateDoubleMatrix(MBSystemClass::nStates, MBSystemClass::nStates, mxREAL);
        for(int irow=0; irow<MBSystemClass::nStates; ++irow)
            for(int icol=0; icol<MBSystemClass::nStates; ++icol)
                mxGetPr(plhs[3])[irow + icol*MBSystemClass::nStates]= system.K(irow, icol);
    }
}
