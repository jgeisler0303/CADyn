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

bool tryGetOption(double &value, const char *name, const mxArray *mxOptions);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    if(nrhs<4 || nrhs>6) { mexErrMsgIdAndTxt("CADyn:InvalidArgument", "Wrong number of arguments. Expecting (x0, dx0, u, p, {ts, {options}})"); return; }
    if(nlhs<1 || nlhs>8) { mexErrMsgIdAndTxt("CADyn:InvalidArgument", "Wrong number of return values. Expecting [x, {dx, {ddx, {converged, {cpu_time, {error, {n_steps, {n_back_steps}}}}}}}]"); return; }
    
    if(!mxIsDouble(prhs[0]) || mxGetNumberOfElements(prhs[0])!=MBSystemClass::nStates) { mexErrMsgIdAndTxt("CADyn:InvalidArgument", "Wrong number of elements in 'x0' (%d expected)", MBSystemClass::nStates); return; }
    if(!mxIsDouble(prhs[1]) || mxGetNumberOfElements(prhs[1])!=MBSystemClass::nStates) { mexErrMsgIdAndTxt("CADyn:InvalidArgument", "Wrong number of elements in 'dx0' (%d expected)", MBSystemClass::nStates); return; }
    // TODO: enable more than on sim step
    if(!mxIsDouble(prhs[2]) || mxGetNumberOfElements(prhs[1])!=MBSystemClass::nInputs) { mexErrMsgIdAndTxt("CADyn:InvalidArgument", "Wrong number of elements in 'u' (%d expected)", MBSystemClass::nInputs); return; }
    
    if(nrhs>=5)
        if(!mxIsDouble(prhs[4]) || mxGetNumberOfElements(prhs[3])!=1) { mexErrMsgIdAndTxt("CADyn:InvalidArgument", "Wrong number of elements in 'ts' (1 expected)"); return; }
    
    const mxArray *mxParams= prhs[3];
    if(!mxIsStruct(mxParams) || mxGetNumberOfElements(mxParams)!=1) {
        mexErrMsgIdAndTxt("CADyn:InvalidArgument", "Input p must be a scalar struct.\n");
        return;
    }
    
    if(nrhs>=6) {
        const mxArray *mxOptions= prhs[5];
        if(!mxIsStruct(mxOptions) || mxGetNumberOfElements(mxOptions)!=1) {
            mexErrMsgIdAndTxt("CADyn:InvalidArgument", "Input options must be a scalar struct.\n");
            return;
        }
    }
    
    MBSystemClass system;
    
    for(const auto &iter : system.param.info_map) {
        const mxArray *mxParam;
        if((mxParam= mxGetField(mxParams, 0, iter.first.c_str()))==NULL) {
            mexErrMsgIdAndTxt("CADyn:InvalidArgument", "Required parameter '%s' is not member of parameters struct.\n", iter.first.c_str());
            return;
        }
        // for now, only scalar parameters supported
        int m_= mxGetM(mxParam);
        int n_= mxGetN(mxParam);
        if(mxIsSparse(mxParam) || !mxIsDouble(mxParam) || (m_!=1 && n_!=1)) {
            mexErrMsgIdAndTxt("CADyn:InvalidArgument", "Parameter name '%s' must be a vector length %d.\n", iter.first.c_str(), 1);
            return;
        }
        system.param.setParam(iter.first.c_str(), mxGetPr(mxParam)[0]);
    }
    
    if(nrhs>=6) {
        const mxArray *mxOptions= prhs[5];
        double value;
        
        if(tryGetOption(value, "AbsTol", mxOptions))
            system.AbsTol= value;
        
        if(tryGetOption(value, "RelTol", mxOptions))
            system.RelTol= value;

        if(tryGetOption(value, "StepTol", mxOptions))
            system.StepTol= value;

        if(tryGetOption(value, "hminmin", mxOptions))
            system.hminmin= value;

        if(tryGetOption(value, "jac_recalc_step", mxOptions))
            system.jac_recalc_step= value;

        if(tryGetOption(value, "max_steps", mxOptions))
            system.max_steps= value;
        //     doflocked
    }
    
    double *x0= mxGetPr(prhs[0]);
    for(int i=0; i<MBSystemClass::nStates; ++i)
        system.q(i)= x0[i];
    
    double *dx0= mxGetPr(prhs[1]);
    for(int i=0; i<MBSystemClass::nStates; ++i)
        system.qd(i)= dx0[i];
    
    double *u_sim= mxGetPr(prhs[2]);;
    for(int i=0; i<MBSystemClass::nInputs; ++i)
        system.u(i)= u_sim[i];
    
    double ts= 0.0;;
    if(nrhs>=5)
        ts= mxGetScalar(prhs[4]);
    
    std::clock_t startcputime = std::clock();
    
    bool res;
    if(ts==0.0)
        res= system.staticEquilibrium();
    else
        res= system.newmarkIntegration(ts, ts, ts, nullptr);
    
    double cpu_duration = (std::clock() - startcputime) / (double)CLOCKS_PER_SEC;
    VecX q, qd, qdd;
    
    plhs[0]= mxCreateDoubleMatrix(MBSystemClass::nStates, 1, mxREAL);
    for(int i=0; i<MBSystemClass::nStates; ++i)
        mxGetPr(plhs[0])[i]= system.q(i);
        
    if(nlhs>1) {
        plhs[1]= mxCreateDoubleMatrix(MBSystemClass::nStates, 1, mxREAL);
        for(int i=0; i<MBSystemClass::nStates; ++i)
            mxGetPr(plhs[1])[i]= system.qd(i);
    }
    
    if(nlhs>2) {
        plhs[2]= mxCreateDoubleMatrix(MBSystemClass::nStates, 1, mxREAL);
        for(int i=0; i<MBSystemClass::nStates; ++i)
            mxGetPr(plhs[2])[i]= system.qdd(i);
    }
    
    if(nlhs>3) {
        plhs[3]= mxCreateDoubleMatrix(1, 1, mxREAL);
        mxGetPr(plhs[3])[0]= res;
    } else if(!res)
        mexErrMsgIdAndTxt("CADyn:Integrator", "Error in integrator");
        
    if(nlhs>4) {
        plhs[4]= mxCreateDoubleMatrix(1, 1, mxREAL);
        mxGetPr(plhs[4])[0]= cpu_duration;
    }
    
    if(nlhs>5) {
        plhs[5]= mxCreateDoubleMatrix(1, 1, mxREAL);
        mxGetPr(plhs[5])[0]= system.errq;
    }
    
    if(nlhs>6) {
        plhs[6]= mxCreateDoubleMatrix(1, 1, mxREAL);
        mxGetPr(plhs[6])[0]= system.n_steps;
    }

    if(nlhs>7) {
        plhs[7]= mxCreateDoubleMatrix(1, 1, mxREAL);
        mxGetPr(plhs[7])[0]= system.n_back_steps;
    }
}

bool tryGetOption(double &value, const char *name, const mxArray *mxOptions) {
    const mxArray *mxOption;
    if((mxOption= mxGetField(mxOptions, 0, name))!=NULL) {
        int m_= mxGetM(mxOption);
        int n_= mxGetN(mxOption);
        if(mxIsSparse(mxOption) || !mxIsDouble(mxOption) || (m_!=1 && n_!=1)) {
            mexErrMsgIdAndTxt("CADyn:InvalidArgument", "Option name '%s' must be a scalar.\n", name);
            return false;
        }
        value= mxGetPr(mxOption)[0];
        
        return true;
    }
    return false;
}
