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
#include "RK1condensed.hpp"
#include "OptionInfo.hpp"

#include <Eigen/Dense>
#include <Eigen/Geometry>
//#include "T1_est_ode1.hpp" // uncomment to make IntelliSense happy

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

#include "EKF_RK1_autotune.hpp"

enum input_idx {
    in_idx_x0= 0,
    in_idx_u,
    in_idx_y,
    in_idx_param,
    in_idx_ts,
    in_idx_x_ul,
    in_idx_x_ll,
    in_idx_Q,
    in_idx_R,
    in_idx_T_adapt,
    in_idx_adaptScale,
    in_idx_fixedQxx,
    in_idx_fixedRxx,
    in_idx_options,
    in_idx_P,
    in_idx_dx0,
    in_idx_last
};

enum output_idx {
    out_idx_x= 0,
    out_idx_dx,
    out_idx_y,
    out_idx_Q,
    out_idx_R,
    out_idx_time,
    out_idx_d_norm,
    out_idx_p_xx,
    out_idx_r_xx,
    out_idx_s_xx,
    out_idx_P,
    out_idx_last
};

typedef EKF_RK1_autotune<ProblemDefinition> EKF;

bool tryGetOption(double *value, const char *name, const mxArray *mxOptions, int m__=1, int n__=1);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    if(nrhs<(in_idx_x_ll+1) || nrhs>in_idx_last) { mexErrMsgIdAndTxt("CADyn:InvalidArgument", "Wrong number of arguments. Expecting (x0, u, y, param, ts, x_ul, x_ll, {Q, R, T_adapt, adaptScale, fixedQxx, fixedRxx, options, P0, dx0})"); return; }
    if(nlhs<(out_idx_y+1) || nlhs>out_idx_last) { mexErrMsgIdAndTxt("CADyn:InvalidArgument", "Wrong number of return values. Expecting [x, xd, y, {Q, R, cpu_time, d_norm, p_xx, r_xx, s_xx, P_end}]"); return; }
    
    if(!mxIsDouble(prhs[in_idx_x0]) || mxGetNumberOfElements(prhs[in_idx_x0])!=EKF::dimX) { mexErrMsgIdAndTxt("CADyn:InvalidArgument", "Wrong number of elements in 'x0' (%d expected)", EKF::dimX); return; }
    
    if(!mxIsDouble(prhs[in_idx_u]) || mxGetM(prhs[in_idx_u])!=EKF::dimIn) { mexErrMsgIdAndTxt("CADyn:InvalidArgument", "Wrong number of rows in 'u' (%d expected)", EKF::dimIn); return; }
    if(!mxIsDouble(prhs[in_idx_u]) || mxGetN(prhs[in_idx_u])<2) { mexErrMsgIdAndTxt("CADyn:InvalidArgument", "Wrong number of columns in 'u' (at least %d expected)", 2); return; }
    
    if(!mxIsDouble(prhs[in_idx_y]) || mxGetM(prhs[in_idx_y])!=EKF::dimOut) { mexErrMsgIdAndTxt("CADyn:InvalidArgument", "Wrong number of rows in 'y' (%d expected)", EKF::dimOut); return; }
    if(!mxIsDouble(prhs[in_idx_y]) || mxGetN(prhs[in_idx_y])!=mxGetN(prhs[in_idx_u])) { mexErrMsgIdAndTxt("CADyn:InvalidArgument", "Wrong number of columns in 'y' (expected same as in u)"); return; }
    
    if(!mxIsDouble(prhs[in_idx_ts]) || mxGetNumberOfElements(prhs[in_idx_ts])!=1) { mexErrMsgIdAndTxt("CADyn:InvalidArgument", "Wrong number of elements in 'ts' (1 expected)"); return; }
    
    if(!mxIsDouble(prhs[in_idx_x_ul]) || mxGetNumberOfElements(prhs[in_idx_x_ul])!=EKF::dimX) { mexErrMsgIdAndTxt("CADyn:InvalidArgument", "Wrong number of elements in 'x_ul' (%d expected)", EKF::dimX); return; }

    if(!mxIsDouble(prhs[in_idx_x_ll]) || mxGetNumberOfElements(prhs[in_idx_x_ll])!=EKF::dimX) { mexErrMsgIdAndTxt("CADyn:InvalidArgument", "Wrong number of elements in 'x_ll' (%d expected)", EKF::dimX); return; }

    if(nrhs>in_idx_Q) if(!mxIsEmpty(prhs[in_idx_Q])) if(!mxIsDouble(prhs[in_idx_Q]) || mxGetM(prhs[in_idx_Q])!=EKF::dimX || mxGetN(prhs[in_idx_Q])!=EKF::dimX) { mexErrMsgIdAndTxt("CADyn:InvalidArgument", "Wrong number of rows/columns in 'Q' (%d/%d expected)", EKF::dimX, EKF::dimX); return; }
    
    if(nrhs>in_idx_R) if(!mxIsEmpty(prhs[in_idx_R])) if(!mxIsDouble(prhs[in_idx_R]) || mxGetM(prhs[in_idx_R])!=EKF::dimOut || mxGetN(prhs[in_idx_R])!=EKF::dimOut) { mexErrMsgIdAndTxt("CADyn:InvalidArgument", "Wrong number of rows/columns in 'R' (%d/%d expected)", EKF::dimOut, EKF::dimOut); return; }
    
    if(nrhs>in_idx_T_adapt) if(!mxIsDouble(prhs[in_idx_T_adapt]) || mxGetM(prhs[in_idx_T_adapt])!=1 || mxGetN(prhs[in_idx_T_adapt])!=1) { mexErrMsgIdAndTxt("CADyn:InvalidArgument", "Wrong number of elements in 'in_idx_T_adapt' (1 expected)"); return; }

    if(!mxIsDouble(prhs[in_idx_adaptScale]) || mxGetNumberOfElements(prhs[in_idx_adaptScale])!=EKF::dimOut) { mexErrMsgIdAndTxt("CADyn:InvalidArgument", "Wrong number of elements in 'in_idx_adaptScale' (%d expected)", EKF::dimOut); return; }
    
    if(!mxIsDouble(prhs[in_idx_fixedQxx]) || mxGetNumberOfElements(prhs[in_idx_fixedQxx])!=EKF::dimX) { mexErrMsgIdAndTxt("CADyn:InvalidArgument", "Wrong number of elements in 'in_idx_fixedQxx' (%d expected)", EKF::dimX); return; }

    if(!mxIsDouble(prhs[in_idx_fixedRxx]) || mxGetNumberOfElements(prhs[in_idx_fixedRxx])!=EKF::dimOut) { mexErrMsgIdAndTxt("CADyn:InvalidArgument", "Wrong number of elements in 'in_idx_fixedRxx' (%d expected)", EKF::dimOut); return; }

    const mxArray *mxParams= prhs[in_idx_param];
    if(!mxIsStruct(mxParams) || mxGetNumberOfElements(mxParams)!=1) {
        mexErrMsgIdAndTxt("CADyn:InvalidArgument", "Input param must be a scalar struct.\n");
        return;
    }
    
    if(nrhs>in_idx_options) {
        const mxArray *mxOptions= prhs[in_idx_options];
        if(!mxIsStruct(mxOptions) || mxGetNumberOfElements(mxOptions)!=1) {
            mexErrMsgIdAndTxt("CADyn:InvalidArgument", "Input options must be a scalar struct.\n");
            return;
        }
    }
    if(nrhs>in_idx_dx0)
        if(!mxIsDouble(prhs[in_idx_dx0]) || mxGetNumberOfElements(prhs[in_idx_dx0])!=EKF::dimX) { mexErrMsgIdAndTxt("CADyn:InvalidArgument", "Wrong number of elements in 'dx0' (%d expected)", EKF::dimX); return; }

    EKF ekf;

    for(const auto &iter : ekf.rk1_solver.param.info_map) {
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
        ekf.rk1_solver.param.setParam(iter.first.c_str(), mxGetPr(mxParam));
    }
    
    if(nrhs>in_idx_options) {
        const mxArray *mxOptions= prhs[in_idx_options];
        try {
            for(const auto &iter : ekf.rk1_solver.options_map) {
                double value;
                if(tryGetOption(&value, iter.first.c_str(), mxOptions)) {
                    ekf.rk1_solver.setOpt(iter.first.c_str(), value);
                }
            }
        } catch (const std::exception& e) {
            mexWarnMsgIdAndTxt("CADyn:Integrator", "CADyn Error setting option: %s", e.what());
        }

        double *locked_states= (double*)mxMalloc(ekf.rk1_solver.nX*sizeof(double));
        if(tryGetOption(locked_states, "locked_states", mxOptions, ekf.rk1_solver.nX, 1)) {
            for(int i= 0; i<ekf.rk1_solver.nX; i++)
                ekf.rk1_solver.locked_states[i]= locked_states[i]!=0.0;
        }
        mxFree(locked_states);
    }
    
    double *x0= mxGetPr(prhs[in_idx_x0]);
    for(int i= 0; i<EKF::dimX; ++i) {
        ekf.rk1_solver.x(i)= x0[i];
    }
    if(nrhs>in_idx_dx0) {
        double *dx0= mxGetPr(prhs[in_idx_dx0]);
        for(int i=ekf.dimI; i<ekf.rk1_solver.nX; ++i)
            ekf.rk1_solver.xdot(i-ekf.dimI)= dx0[i];
            // TODO: check if dx0.head(nI) == E*x0.tail(nII)
    }
    
    double *x_ul= mxGetPr(prhs[in_idx_x_ul]);
    for(int i=0; i<EKF::dimX; ++i)
        ekf.x_ul(i)= x_ul[i];
    double *x_ll= mxGetPr(prhs[in_idx_x_ll]);
    for(int i=0; i<EKF::dimX; ++i)
        ekf.x_ll(i)= x_ll[i];

    if(nrhs>in_idx_Q && !mxIsEmpty(prhs[in_idx_Q])) {
        double *Q= mxGetPr(prhs[in_idx_Q]);
        for(int i=0; i<EKF::dimX; ++i)
            for(int j=0; j<EKF::dimX; ++j)
                ekf.ekfQ(i, j)= Q[i + j*EKF::dimX];
    }
    if(nrhs>in_idx_R && !mxIsEmpty(prhs[in_idx_R])) {
        double *R= mxGetPr(prhs[in_idx_R]);
        for(int i=0; i<EKF::dimOut; ++i)
            for(int j=0; j<EKF::dimOut; ++j)
                ekf.ekfR(i, j)= R[i + j*EKF::dimOut];
    }
    if(nrhs>in_idx_P && !mxIsEmpty(prhs[in_idx_P])) {
        double *P= mxGetPr(prhs[in_idx_P]);
        for(int i=0; i<EKF::dimX; ++i)
            for(int j=0; j<EKF::dimX; ++j)
                ekf.ekfSigma(i, j)= P[i + j*EKF::dimX];
    }
    if(nrhs>in_idx_T_adapt)
        ekf.T_adapt= mxGetScalar(prhs[in_idx_T_adapt]);
        
    if(nrhs>in_idx_adaptScale && !mxIsEmpty(prhs[in_idx_adaptScale])) {
        double *adaptScale= mxGetPr(prhs[in_idx_adaptScale]);
        for(int i=0; i<EKF::dimOut; ++i)
            ekf.adaptScale(i)= adaptScale[i];
    }

    if(nrhs>in_idx_fixedQxx && !mxIsEmpty(prhs[in_idx_fixedQxx])) {
        double *fixedQxx= mxGetPr(prhs[in_idx_fixedQxx]);
        for(int i=0; i<EKF::dimX; ++i)
            ekf.fixedQxx(i)= fixedQxx[i];
    }

    if(nrhs>in_idx_fixedRxx && !mxIsEmpty(prhs[in_idx_fixedRxx])) {
        double *fixedRxx= mxGetPr(prhs[in_idx_fixedRxx]);
        for(int i=0; i<EKF::dimOut; ++i)
            ekf.fixedRxx(i)= fixedRxx[i];
    }

    double *u= mxGetPr(prhs[in_idx_u]);
    double *y_meas= mxGetPr(prhs[in_idx_y]);
    
    double ts= mxGetScalar(prhs[in_idx_ts]);
    
    ekf.rk1_solver.precalcConsts();

    plhs[out_idx_x]= mxCreateDoubleMatrix(EKF::dimX, mxGetN(prhs[in_idx_u]), mxREAL);
    plhs[out_idx_dx]= mxCreateDoubleMatrix(EKF::dimX, mxGetN(prhs[in_idx_u]), mxREAL);
    plhs[out_idx_y]= mxCreateDoubleMatrix(EKF::dimOut, mxGetN(prhs[in_idx_u]), mxREAL);

    if(nlhs>out_idx_d_norm) {
        plhs[out_idx_d_norm]= mxCreateDoubleMatrix(1, mxGetN(prhs[in_idx_u]), mxREAL);
    }    
    if(nlhs>out_idx_p_xx) {
        plhs[out_idx_p_xx]= mxCreateDoubleMatrix(EKF::dimX, mxGetN(prhs[in_idx_u]), mxREAL);
    }    
    if(nlhs>out_idx_r_xx) {
        plhs[out_idx_r_xx]= mxCreateDoubleMatrix(EKF::dimOut, mxGetN(prhs[in_idx_u]), mxREAL);
    }    
    if(nlhs>out_idx_s_xx) {
        plhs[out_idx_s_xx]= mxCreateDoubleMatrix(EKF::dimOut, mxGetN(prhs[in_idx_u]), mxREAL);
    }    
    
    EKF::VecOut y_meas_vec;
    
    std::clock_t startcputime = std::clock();
    try {
        for(int i= 0; i<(int)mxGetN(prhs[in_idx_u]); ++i) {
            if(i>0) {
                for(int j= 0; j<ekf.dimIn; ++j)
                    ekf.rk1_solver.u(j)= u[j + (i-1)*ekf.dimIn];
                
                for(int j= 0; j<ekf.dimOut; ++j)
                    y_meas_vec(j)= y_meas[j + i*ekf.dimOut];
                    
                ekf.next(ts, y_meas_vec);
            }
            for(int j=0; j<EKF::dimX; ++j)
                mxGetPr(plhs[out_idx_x])[j + i*EKF::dimX]= ekf.rk1_solver.x(j);
            
            ProblemDefinition::VecI dxI = ekf.rk1_solver.E * ekf.rk1_solver.x.tail(ekf.rk1_solver.nII);
            for(int j=0; j<EKF::dimI; ++j)
                mxGetPr(plhs[out_idx_dx])[j + i*EKF::dimX]= dxI(j);
            for(int j=EKF::dimI; j<EKF::dimX; ++j)
                mxGetPr(plhs[out_idx_dx])[j + i*EKF::dimX]= ekf.rk1_solver.xdot(j-EKF::dimI);
            for(int j=0; j<EKF::dimOut; ++j)
                mxGetPr(plhs[out_idx_y])[j + i*EKF::dimOut]= ekf.rk1_solver.y(j);
            
            if(nlhs>out_idx_d_norm) {
                mxGetPr(plhs[out_idx_d_norm])[i]= ekf.d_norm;
            }    
            if(nlhs>out_idx_p_xx) {
                for(int j=0; j<EKF::dimX; ++j)
                    mxGetPr(plhs[out_idx_p_xx])[j + i*EKF::dimX]= ekf.ekfSigma_pred(j, j);
            }    
            if(nlhs>out_idx_r_xx) {
                for(int j=0; j<EKF::dimOut; ++j)
                    mxGetPr(plhs[out_idx_r_xx])[j + i*EKF::dimOut]= ekf.ekfR(j, j);
            }    
            if(nlhs>out_idx_s_xx) {
                for(int j=0; j<EKF::dimOut; ++j)
                mxGetPr(plhs[out_idx_s_xx])[j + i*EKF::dimOut]= ekf.ekfS(j, j);
            }    
        }
    } catch(std::exception &e) {
        mexWarnMsgIdAndTxt("CADyn:EKF", "Error EKF: %s", e.what());
    }
    double cpu_duration = (std::clock() - startcputime) / (double)CLOCKS_PER_SEC;
    
        
    if(nlhs>out_idx_Q) {
        plhs[out_idx_Q]= mxCreateDoubleMatrix(EKF::dimX, EKF::dimX, mxREAL);
        for(int i=0; i<EKF::dimX; ++i)
            for(int j=0; j<EKF::dimX; ++j)
                mxGetPr(plhs[out_idx_Q])[j + i*EKF::dimX]= ekf.ekfQ(j, i);
    }    
    if(nlhs>out_idx_R) {
        plhs[out_idx_R]= mxCreateDoubleMatrix(EKF::dimOut, EKF::dimOut, mxREAL);
        for(int i=0; i<EKF::dimOut; ++i)
            for(int j=0; j<EKF::dimOut; ++j)
                mxGetPr(plhs[out_idx_R])[i + j*EKF::dimOut]= ekf.ekfR(i, j);
    }    
    if(nlhs>out_idx_P) {
        plhs[out_idx_P]= mxCreateDoubleMatrix(EKF::dimX, EKF::dimX, mxREAL);
        for(int i=0; i<EKF::dimX; ++i)
            for(int j=0; j<EKF::dimX; ++j)
                mxGetPr(plhs[out_idx_P])[i + j*EKF::dimX]= ekf.ekfSigma(i, j);
    }    
    if(nlhs>out_idx_time) {
        plhs[out_idx_time]= mxCreateDoubleMatrix(1, 1, mxREAL);
        mxGetPr(plhs[out_idx_time])[0]= cpu_duration;
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
