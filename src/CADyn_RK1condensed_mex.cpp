#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <exception>
#include <iosfwd>
#include <cmath>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <algorithm>
#include <ctime>
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include "OptionInfo.hpp"
#include "RK1condensed.hpp"

#define makeProblemDefinition(c) c
#define exp_makeProblemDefinition(c) makeProblemDefinition(c)
#define ProblemDefinition exp_makeProblemDefinition(MBSystem)

#define stringify(x)  #x
#define expand_and_stringify(x) stringify(x ## _ode1.hpp)
#define INCL_FILE_STR(x) expand_and_stringify(x)
#include INCL_FILE_STR(MBSystem)

#include "mex.h"
#ifndef  HAVE_OCTAVE
#include "matrix.h"
#endif

enum input_idx {
    in_idx_x0= 0,
    in_idx_dx0,
    in_idx_u,
    in_idx_param,
    in_idx_ts,
    in_idx_options,
    in_idx_last
};

enum output_idx {
    out_idx_x= 0,
    out_idx_dx,
    out_idx_y,
    out_idx_fxfu,
    out_idx_gxgu,
    out_idx_converged,
    out_idx_cpu_time,
    out_idx_error,
    out_idx_n_steps,
    out_idx_n_back_steps,
    out_idx_sub_iter,
    out_idx_last
};

bool tryGetOption(double *value, const char *name, const mxArray *mxOptions, int m__=1, int n__=1);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    if(nrhs<in_idx_ts || nrhs>in_idx_last) { mexErrMsgIdAndTxt("CADyn:InvalidArgument", "Wrong number of arguments. Expecting (x0, dx0, u, p, {ts, {options}})"); return; }
    if(nlhs<1 || nlhs>out_idx_last) { mexErrMsgIdAndTxt("CADyn:InvalidArgument", "Wrong number of return values. Expecting [x, {dx, {y, {fxfu, {gxgu, {converged, {cpu_time, {error, {n_steps, {n_back_steps, {sub_iter}}}}}}}}}}}]"); return; }
    if(nrhs==4 && nlhs!=1) { mexErrMsgIdAndTxt("CADyn:InvalidArgument", "When four arguments are supplied only the output is calculated and returned as the one and only return value"); return; }
    
    if(!mxIsDouble(prhs[in_idx_x0]) || mxGetNumberOfElements(prhs[in_idx_x0])!=ProblemDefinition::dimX) { mexErrMsgIdAndTxt("CADyn:InvalidArgument", "Wrong number of elements in 'x0' (%d expected)", ProblemDefinition::dimX); return; }
    if(!mxIsDouble(prhs[in_idx_dx0]) || mxGetNumberOfElements(prhs[in_idx_dx0])!=ProblemDefinition::dimX) { mexErrMsgIdAndTxt("CADyn:InvalidArgument", "Wrong number of elements in 'dx0' (%d expected)", ProblemDefinition::dimX); return; }
    if(!mxIsDouble(prhs[in_idx_u]) || mxGetNumberOfElements(prhs[in_idx_u])!=ProblemDefinition::dimIn) { mexErrMsgIdAndTxt("CADyn:InvalidArgument", "Wrong number of elements in 'u' (%d expected)", ProblemDefinition::dimIn); return; }
    
    const mxArray *mxParams= prhs[in_idx_param];
    if(!mxIsStruct(mxParams) || mxGetNumberOfElements(mxParams)!=1) {
        mexErrMsgIdAndTxt("CADyn:InvalidArgument", "Input p must be a scalar struct.\n");
        return;
    }
    
    if(nrhs>in_idx_ts)
        if(!mxIsDouble(prhs[in_idx_ts]) || mxGetNumberOfElements(prhs[in_idx_ts])!=1) { mexErrMsgIdAndTxt("CADyn:InvalidArgument", "Wrong number of elements in 'ts' (1 expected)"); return; }
    
    if(nrhs>in_idx_options) {
        const mxArray *mxOptions= prhs[in_idx_options];
        if(!mxIsStruct(mxOptions) || mxGetNumberOfElements(mxOptions)!=1) {
            mexErrMsgIdAndTxt("CADyn:InvalidArgument", "Input options must be a scalar struct.\n");
            return;
        }
    }
    
    RK1condensed<ProblemDefinition> rk1;
    
    for(const auto &iter : rk1.param.info_map) {
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
        rk1.param.setParam(iter.first.c_str(), mxGetPr(mxParam));
    }
    rk1.precalcConsts();
    
    if(nrhs>in_idx_options) {
        const mxArray *mxOptions= prhs[in_idx_options];
        try {
            for(const auto &iter : rk1.options_map) {
                double value;
                if(tryGetOption(&value, iter.first.c_str(), mxOptions)) {
                    rk1.setOpt(iter.first.c_str(), value);
                }
            }
        } catch (const std::exception& e) {
            mexWarnMsgIdAndTxt("CADyn:Integrator", "CADyn Error setting option: %s", e.what());
        }

        double *locked_states= (double*)mxMalloc(rk1.nX*sizeof(double));
        if(tryGetOption(locked_states, "locked_states", mxOptions, rk1.nX, 1)) {
            for(int i= 0; i<rk1.nX; i++)
                rk1.locked_states[i]= locked_states[i]!=0.0;
        }
        mxFree(locked_states);
    }
    
    double *x0= mxGetPr(prhs[in_idx_x0]);
    for(int i=0; i<rk1.nX; ++i)
        rk1.x(i)= x0[i];
    
    double *dx0= mxGetPr(prhs[in_idx_dx0]);
    for(int i=rk1.nI; i<rk1.nX; ++i)
        rk1.xdot(i-rk1.nI)= dx0[i];
        // TODO: check if dx0.head(nI) == E*x0.tail(nII)
    
    double *u_sim= mxGetPr(prhs[in_idx_u]);;
    for(int i=0; i<rk1.nIn; ++i)
        rk1.u(i)= u_sim[i];
    
    double ts= 0.0;
    if(nrhs>in_idx_ts) {
        ts= mxGetScalar(prhs[in_idx_ts]);
    } else {
        // No ts means: get outputs
        rk1.evaluateOutput();
        plhs[0]= mxCreateDoubleMatrix(rk1.nOut, 1, mxREAL);
        for(int i=0; i<rk1.nOut; ++i)
            mxGetPr(plhs[0])[i]= rk1.y(i);
        
        return;
    }
    
    std::clock_t startcputime = std::clock();
    
    bool res= true;
    try {
        if(ts==-1)
            rk1.staticEquilibriumWithLin();
        else
            rk1.intervalWithSens(ts);
    } catch (const std::exception& e) {
        mexWarnMsgIdAndTxt("CADyn:Integrator", "CADyn Error: %s", e.what());
        res= false;
    }
    
    double cpu_duration = (std::clock() - startcputime) / (double)CLOCKS_PER_SEC;
    
    plhs[out_idx_x]= mxCreateDoubleMatrix(rk1.nX, 1, mxREAL);
    for(int i=0; i<rk1.nX; ++i)
        mxGetPr(plhs[out_idx_x])[i]= rk1.x(i);
        
    if(nlhs>out_idx_dx) {
        plhs[out_idx_dx]= mxCreateDoubleMatrix(rk1.nX, 1, mxREAL);
        for(int i=0; i<rk1.nI; ++i)
            mxGetPr(plhs[out_idx_dx])[i]= rk1.E.row(i) * rk1.x.tail(rk1.nII);
            
        for(int i=rk1.nI; i<rk1.nX; ++i)
            mxGetPr(plhs[out_idx_dx])[i]= rk1.xdot(i-rk1.nI);
    }
    
    if(nlhs>out_idx_y) {
        rk1.evaluateOutput();
        plhs[out_idx_y]= mxCreateDoubleMatrix(rk1.nOut, 1, mxREAL);
        for(int i=0; i<rk1.nOut; ++i)
            mxGetPr(plhs[out_idx_y])[i]= rk1.y(i);
    }
    
    if(nlhs>out_idx_fxfu) {
        plhs[out_idx_fxfu]= mxCreateDoubleMatrix(rk1.nX, rk1.nX+rk1.nIn, mxREAL);
        for(int irow=0; irow<rk1.nX; ++irow) {
            for(int icol=0; icol<rk1.nX; ++icol)
                mxGetPr(plhs[out_idx_fxfu])[irow + icol*rk1.nX]= rk1.fx(irow, icol);
            for(int icol=0; icol<rk1.nIn; ++icol)
                mxGetPr(plhs[out_idx_fxfu])[irow + (icol+rk1.nX)*rk1.nX]= rk1.fu(irow, icol);
        }
    }
    
    if(nlhs>out_idx_gxgu) {
        plhs[out_idx_gxgu]= mxCreateDoubleMatrix(rk1.nOut, rk1.nX+rk1.nIn, mxREAL);
        for(int irow=0; irow<rk1.nOut; ++irow) {
            for(int icol=0; icol<rk1.nX; ++icol)
                mxGetPr(plhs[out_idx_gxgu])[irow + icol*rk1.nOut]= rk1.gx(irow, icol);
            for(int icol=0; icol<rk1.nIn; ++icol)
                mxGetPr(plhs[out_idx_gxgu])[irow + (icol+rk1.nX)*rk1.nOut]= rk1.gu(irow, icol);
        }
    }
    
    if(nlhs>out_idx_converged) {
        plhs[out_idx_converged]= mxCreateDoubleMatrix(1, 1, mxREAL);
        mxGetPr(plhs[out_idx_converged])[0]= res;
    }
        
    if(nlhs>out_idx_cpu_time) {
        plhs[out_idx_cpu_time]= mxCreateDoubleMatrix(1, 1, mxREAL);
        mxGetPr(plhs[out_idx_cpu_time])[0]= cpu_duration;
    }
    
    if(nlhs>out_idx_error) {
        plhs[out_idx_error]= mxCreateDoubleMatrix(1, 1, mxREAL);
        mxGetPr(plhs[out_idx_error])[0]= rk1.computeError(ts);
    }
    
    if(nlhs>out_idx_n_steps) {
        plhs[out_idx_n_steps]= mxCreateDoubleMatrix(1, 1, mxREAL);
        mxGetPr(plhs[out_idx_n_steps])[0]= rk1.n_steps;
    }

    if(nlhs>out_idx_n_back_steps) {
        plhs[out_idx_n_back_steps]= mxCreateDoubleMatrix(1, 1, mxREAL);
        mxGetPr(plhs[out_idx_n_back_steps])[0]= rk1.n_back_steps;
    }   
    
    if(nlhs>out_idx_sub_iter) {
        plhs[out_idx_sub_iter]= mxCreateDoubleMatrix(1, 1, mxREAL);
        mxGetPr(plhs[out_idx_sub_iter])[0]= rk1.sub_iter;
    }   
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
