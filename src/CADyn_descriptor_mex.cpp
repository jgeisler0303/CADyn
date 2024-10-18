#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define makeMBSystemClass(c) c
#define exp_makeMBSystemClass(c) makeMBSystemClass(c)
#define MBSystemClass exp_makeMBSystemClass(MBSystem)

#define stringify(x)  #x
#define expand_and_stringify(x) stringify(x ## _descriptor_form.hpp)
#define INCL_FILE_STR(x) expand_and_stringify(x)

#include INCL_FILE_STR(MBSystem)

#include "mex.h"
#ifndef  HAVE_OCTAVE
#include "matrix.h"
#endif

enum input_idx {
    in_idx_q0= 0,
    in_idx_dq0,
    in_idx_u,
    in_idx_param,
    in_idx_last
};

enum output_idx {
    out_idx_M= 0,
    out_idx_f,
    out_idx_A,
    out_idx_B,    
    out_idx_C,
    out_idx_D,
    out_idx_last
};

bool tryGetOption(double *value, const char *name, const mxArray *mxOptions, int m__=1, int n__=1);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    if(nrhs!=in_idx_last) { mexErrMsgIdAndTxt("CADyn:InvalidArgument", "Wrong number of arguments. Expecting (q0, dq0, u, p)"); return; }
    if(nlhs!=out_idx_last) { mexErrMsgIdAndTxt("CADyn:InvalidArgument", "Wrong number of return values. Expecting [M, f, A, B, C, D]"); return; }
    
    if(!mxIsDouble(prhs[in_idx_q0]) || mxGetNumberOfElements(prhs[in_idx_q0])!=MBSystemClass::nbrdof) { mexErrMsgIdAndTxt("CADyn:InvalidArgument", "Wrong number of elements in 'q0' (%d expected)", MBSystemClass::nbrdof); return; }
    
    if(!mxIsDouble(prhs[in_idx_dq0]) || mxGetNumberOfElements(prhs[in_idx_dq0])!=MBSystemClass::nbrdof) { mexErrMsgIdAndTxt("CADyn:InvalidArgument", "Wrong number of elements in 'dd0' (%d expected)", MBSystemClass::nbrdof); return; }
    
    if(!mxIsDouble(prhs[in_idx_u]) || mxGetNumberOfElements(prhs[in_idx_u])!=MBSystemClass::nbrin) { mexErrMsgIdAndTxt("CADyn:InvalidArgument", "Wrong number of elements in 'u' (%d expected)", MBSystemClass::nbrin); return; }
    
    const mxArray *mxParams= prhs[in_idx_param];
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
    
    double *q0= mxGetPr(prhs[in_idx_q0]);
    for(int i=0; i<MBSystemClass::nbrdof; ++i)
        system.q(i)= q0[i];
    
    double *dq0= mxGetPr(prhs[in_idx_dq0]);
    for(int i=0; i<MBSystemClass::nbrdof; ++i)
        system.qd(i)= dq0[i];
    
    double *u_sim= mxGetPr(prhs[in_idx_u]);;
    for(int i=0; i<MBSystemClass::nbrin; ++i)
        system.u(i)= u_sim[i];
    
    try {
        system.calcMf();
        system.calcJacobian();
        system.calcB();
        system.calcCDF();
        // system.calcOut();
    } catch (const std::exception& e) {
        mexErrMsgIdAndTxt("CADyn:DescriptorForm", "CADyn Error: %s", e.what());
        return;
    }
    
    plhs[out_idx_M]= mxCreateDoubleMatrix(MBSystemClass::nbrdof, MBSystemClass::nbrdof, mxREAL);
    for(int i=0; i<MBSystemClass::nbrdof; ++i)
        for(int j=0; j<MBSystemClass::nbrdof; ++j)
            mxGetPr(plhs[out_idx_M])[i + j*MBSystemClass::nbrdof]= system.M(i, j);
        
    plhs[out_idx_f]= mxCreateDoubleMatrix(MBSystemClass::nbrdof, 1, mxREAL);
    for(int i=0; i<MBSystemClass::nbrdof; ++i)
        mxGetPr(plhs[out_idx_f])[i]= system.f(i);

    plhs[out_idx_A]= mxCreateDoubleMatrix(2*MBSystemClass::nbrdof, 2*MBSystemClass::nbrdof, mxREAL);
    for(int i=0; i<MBSystemClass::nbrdof; ++i)
        for(int j=0; j<MBSystemClass::nbrdof; ++j)
            mxGetPr(plhs[out_idx_A])[i + j*2*MBSystemClass::nbrdof]= 0.0;
    for(int i=0; i<MBSystemClass::nbrdof; ++i)
        for(int j=0; j<MBSystemClass::nbrdof; ++j)
            if(i==j)
                mxGetPr(plhs[out_idx_A])[i + (j+MBSystemClass::nbrdof)*2*MBSystemClass::nbrdof]= 1.0;
            else
                mxGetPr(plhs[out_idx_A])[i + (j+MBSystemClass::nbrdof)*2*MBSystemClass::nbrdof]= 0.0;
    for(int i=0; i<MBSystemClass::nbrdof; ++i)
        for(int j=0; j<MBSystemClass::nbrdof; ++j)
            mxGetPr(plhs[out_idx_A])[i+MBSystemClass::nbrdof + j*2*MBSystemClass::nbrdof]= system.K(i, j);
    for(int i=0; i<MBSystemClass::nbrdof; ++i)
        for(int j=0; j<MBSystemClass::nbrdof; ++j)
            mxGetPr(plhs[out_idx_A])[i+MBSystemClass::nbrdof + (j+MBSystemClass::nbrdof)*2*MBSystemClass::nbrdof]= system.C(i, j);

    plhs[out_idx_B]= mxCreateDoubleMatrix(2*MBSystemClass::nbrdof, MBSystemClass::nbrin, mxREAL);
    for(int i=0; i<MBSystemClass::nbrdof; ++i)
        for(int j=0; j<MBSystemClass::nbrin; ++j)
            mxGetPr(plhs[out_idx_B])[i + j*2*MBSystemClass::nbrdof]= 0.0;
    for(int i= 0; i<MBSystemClass::nbrdof; ++i)
        for(int j=0; j<MBSystemClass::nbrin; ++j)
            mxGetPr(plhs[out_idx_B])[i + MBSystemClass::nbrdof + j*2*MBSystemClass::nbrdof]= system.B(i, j);

    plhs[out_idx_C]= mxCreateDoubleMatrix(MBSystemClass::nbrout, 2*MBSystemClass::nbrdof, mxREAL);
    for(int i=0; i<MBSystemClass::nbrout; ++i)
        for(int j=0; j<2*MBSystemClass::nbrdof; ++j)
            mxGetPr(plhs[out_idx_C])[i + j*MBSystemClass::nbrout]= system.CD(i, j);

    plhs[out_idx_D]= mxCreateDoubleMatrix(MBSystemClass::nbrout, MBSystemClass::nbrin, mxREAL);
    for(int i=0; i<MBSystemClass::nbrout; ++i)
        for(int j=0; j<MBSystemClass::nbrin; ++j)
            mxGetPr(plhs[out_idx_D])[i + j*MBSystemClass::nbrout]= system.CD(i, 2*MBSystemClass::nbrdof + j);
}

bool tryGetOption(double *value, const char *name, const mxArray *mxOptions, int m__, int n__) {
    const mxArray *mxOption;
    if((mxOption= mxGetField(mxOptions, 0, name))!=NULL) {
        int m_= mxGetM(mxOption);
        int n_= mxGetN(mxOption);
        if(mxIsSparse(mxOption) || !mxIsDouble(mxOption) || (m_!=m__ && n_!=n__)) {
            mexErrMsgIdAndTxt("CADyn:InvalidArgument", "Option name '%s' must be a scalar.\n", name);
            return false;
        }
        
        for(int i= 0; i<m__; i++)
            for(int j= 0; j<n__; j++)
                value[i + j*m__]= mxGetPr(mxOption)[i + j*m__];
        
        return true;
    }
    return false;
}
