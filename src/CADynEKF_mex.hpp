#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "mex.h"
#ifndef  HAVE_OCTAVE
#include "matrix.h"
#endif

#include "EKF_autotune.hpp"

#define rhs_idx_q0 0
#define rhs_idx_dq0 1
#define rhs_idx_u 2
#define rhs_idx_y 3
#define rhs_idx_p 4
#define rhs_idx_ts 5
#define rhs_idx_x_ul 6
#define rhs_idx_x_ll 7
#define rhs_idx_Q 8
#define rhs_idx_R 9
#define rhs_idx_T_adapt 10
#define rhs_idx_adaptScale 11
#define rhs_idx_fixedQxx 12
#define rhs_idx_fixedRxx 13
#define rhs_idx_opt 14
#define rhs_idx_P 15

#define lhs_idx_q 0
#define lhs_idx_qd 1
#define lhs_idx_qdd 2
#define lhs_idx_y 3
#define lhs_idx_Q 4
#define lhs_idx_R 5
#define lhs_idx_time 6
#define lhs_idx_d_norm 7
#define lhs_idx_p_xx 8
#define lhs_idx_r_xx 9
#define lhs_idx_s_xx 10
#define lhs_idx_P 11

typedef EKF_autotune<EKF_STATES, SYSTEM> EKF;

bool tryGetOption(double *value, const char *name, const mxArray *mxOptions, int m__=1, int n__=1);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    if(nrhs==0 && nlhs==1) {
        plhs[0]= mxCreateDoubleMatrix(1, EKF::nbrstates, mxREAL);
        {
            int x_idx= 0;
            for(int j= 0; j<(int)(sizeof(estimated_q)/sizeof(estimated_q[0])); ++j) {
                if(x_idx<EKF::nbrstates) {
                    mxGetPr(plhs[0])[x_idx]= estimated_q[j]+1;
                    x_idx++;
                }
            }
            for(int j= 0; j<(int)(sizeof(estimated_dq)/sizeof(estimated_dq[0])); ++j) {
                if(x_idx<EKF::nbrstates) {
                    mxGetPr(plhs[0])[x_idx]= estimated_dq[j] + 1 + EKF::nbrdof;
                    x_idx++;
                }
            }
        }        
        return;
    }
    
    if(nrhs<8 || nrhs>16) { mexErrMsgIdAndTxt("CADyn:InvalidArgument", "Wrong number of arguments. Expecting (q0, dq0, u, y, param, ts, x_ul, x_ll, {Q, R, T_adapt, adaptScale, fixedQxx, fixedRxx, options, P0})"); return; }
    if(nlhs<4 || nlhs>12) { mexErrMsgIdAndTxt("CADyn:InvalidArgument", "Wrong number of return values. Expecting [q, qd, qdd, y, {Q, R, cpu_time, d_norm, p_xx, r_xx, s_xx, P_end}]"); return; }
    
    if(!mxIsDouble(prhs[rhs_idx_q0]) || mxGetNumberOfElements(prhs[rhs_idx_q0])!=EKF::nbrdof) { mexErrMsgIdAndTxt("CADyn:InvalidArgument", "Wrong number of elements in 'q0' (%d expected)", EKF::nbrdof); return; }
    
    if(!mxIsDouble(prhs[rhs_idx_dq0]) || mxGetNumberOfElements(prhs[rhs_idx_dq0])!=EKF::nbrdof) { mexErrMsgIdAndTxt("CADyn:InvalidArgument", "Wrong number of elements in 'dq0' (%d expected)", EKF::nbrdof); return; }
    
    if(!mxIsDouble(prhs[rhs_idx_u]) || mxGetM(prhs[rhs_idx_u])!=EKF::nbrin) { mexErrMsgIdAndTxt("CADyn:InvalidArgument", "Wrong number of rows in 'u' (%d expected)", EKF::nbrin); return; }
    if(!mxIsDouble(prhs[rhs_idx_u]) || mxGetN(prhs[rhs_idx_u])<2) { mexErrMsgIdAndTxt("CADyn:InvalidArgument", "Wrong number of columns in 'u' (at least %d expected)", 2); return; }
    
    if(!mxIsDouble(prhs[rhs_idx_y]) || mxGetM(prhs[rhs_idx_y])!=EKF::nbrout) { mexErrMsgIdAndTxt("CADyn:InvalidArgument", "Wrong number of rows in 'y' (%d expected)", EKF::nbrout); return; }
    if(!mxIsDouble(prhs[rhs_idx_y]) || mxGetN(prhs[rhs_idx_y])!=mxGetN(prhs[rhs_idx_u])) { mexErrMsgIdAndTxt("CADyn:InvalidArgument", "Wrong number of columns in 'y' (expected same as in u)"); return; }
    
    if(!mxIsDouble(prhs[rhs_idx_ts]) || mxGetNumberOfElements(prhs[rhs_idx_ts])!=1) { mexErrMsgIdAndTxt("CADyn:InvalidArgument", "Wrong number of elements in 'ts' (1 expected)"); return; }
    
    if(!mxIsDouble(prhs[rhs_idx_x_ul]) || mxGetNumberOfElements(prhs[rhs_idx_x_ul])!=EKF::nbrstates) { mexErrMsgIdAndTxt("CADyn:InvalidArgument", "Wrong number of elements in 'x_ul' (%d expected)", EKF::nbrstates); return; }

    if(!mxIsDouble(prhs[rhs_idx_x_ll]) || mxGetNumberOfElements(prhs[rhs_idx_x_ll])!=EKF::nbrstates) { mexErrMsgIdAndTxt("CADyn:InvalidArgument", "Wrong number of elements in 'x_ll' (%d expected)", EKF::nbrstates); return; }

    if(nrhs>rhs_idx_Q) if(!mxIsEmpty(prhs[rhs_idx_Q])) if(!mxIsDouble(prhs[rhs_idx_Q]) || mxGetM(prhs[rhs_idx_Q])!=EKF::nbrstates || mxGetN(prhs[rhs_idx_Q])!=EKF::nbrstates) { mexErrMsgIdAndTxt("CADyn:InvalidArgument", "Wrong number of rows/columns in 'Q' (%d/%d expected)", EKF::nbrstates, EKF::nbrstates); return; }
    
    if(nrhs>rhs_idx_R) if(!mxIsEmpty(prhs[rhs_idx_R])) if(!mxIsDouble(prhs[rhs_idx_R]) || mxGetM(prhs[rhs_idx_R])!=EKF::nbrout || mxGetN(prhs[rhs_idx_R])!=EKF::nbrout) { mexErrMsgIdAndTxt("CADyn:InvalidArgument", "Wrong number of rows/columns in 'R' (%d/%d expected)", EKF::nbrout, EKF::nbrout); return; }
    
    if(nrhs>rhs_idx_T_adapt) if(!mxIsDouble(prhs[rhs_idx_T_adapt]) || mxGetM(prhs[rhs_idx_T_adapt])!=1 || mxGetN(prhs[rhs_idx_T_adapt])!=1) { mexErrMsgIdAndTxt("CADyn:InvalidArgument", "Wrong number of elements in 'rhs_idx_T_adapt' (1 expected)"); return; }

    if(!mxIsDouble(prhs[rhs_idx_adaptScale]) || mxGetNumberOfElements(prhs[rhs_idx_adaptScale])!=EKF::nbrout) { mexErrMsgIdAndTxt("CADyn:InvalidArgument", "Wrong number of elements in 'rhs_idx_adaptScale' (%d expected)", EKF::nbrout); return; }
    
    if(!mxIsDouble(prhs[rhs_idx_fixedQxx]) || mxGetNumberOfElements(prhs[rhs_idx_fixedQxx])!=EKF::nbrstates) { mexErrMsgIdAndTxt("CADyn:InvalidArgument", "Wrong number of elements in 'rhs_idx_fixedQxx' (%d expected)", EKF::nbrstates); return; }

    if(!mxIsDouble(prhs[rhs_idx_fixedRxx]) || mxGetNumberOfElements(prhs[rhs_idx_fixedRxx])!=EKF::nbrout) { mexErrMsgIdAndTxt("CADyn:InvalidArgument", "Wrong number of elements in 'rhs_idx_fixedRxx' (%d expected)", EKF::nbrout); return; }

    const mxArray *mxParams= prhs[rhs_idx_p];
    if(!mxIsStruct(mxParams) || mxGetNumberOfElements(mxParams)!=1) {
        mexErrMsgIdAndTxt("CADyn:InvalidArgument", "Input param must be a scalar struct.\n");
        return;
    }
    
    if(nrhs>rhs_idx_opt) {
        const mxArray *mxOptions= prhs[rhs_idx_opt];
        if(!mxIsStruct(mxOptions) || mxGetNumberOfElements(mxOptions)!=1) {
            mexErrMsgIdAndTxt("CADyn:InvalidArgument", "Input options must be a scalar struct.\n");
            return;
        }
    }
    
    EKF ekf;
    {
        int x_idx= 0;
        for(int i= 0; i<EKF::nbrdof; ++i) {
            ekf.qx_idx(i)= -1;
            for(int j= 0; j<(int)(sizeof(estimated_q)/sizeof(estimated_q[0])); ++j) {
                if(i==estimated_q[j]) {
                    ekf.qx_idx(i)= x_idx;
                    x_idx++;
                    break;
                }
            }
        }
        for(int i= 0; i<EKF::nbrdof; ++i) {
            ekf.dqx_idx(i)= -1;
            for(int j= 0; j<(int)(sizeof(estimated_dq)/sizeof(estimated_dq[0])); ++j) {
                if(i==estimated_dq[j]) {
                    ekf.dqx_idx(i)= x_idx;
                    x_idx++;
                    break;
                }
            }
        }
    }    
    
    for(const auto &iter : ekf.system.param.info_map) {
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
        ekf.system.param.setParam(iter.first.c_str(), mxGetPr(mxParam));
    }
    
    if(nrhs>=rhs_idx_opt) {
        const mxArray *mxOptions= prhs[rhs_idx_opt];
        double value;
        
        if(tryGetOption(&value, "AbsTol", mxOptions))
            ekf.system.AbsTol= value;
        
        if(tryGetOption(&value, "RelTol", mxOptions))
            ekf.system.RelTol= value;

        if(tryGetOption(&value, "StepTol", mxOptions))
            ekf.system.StepTol= value;

        if(tryGetOption(&value, "hminmin", mxOptions))
            ekf.system.hminmin= value;

        if(tryGetOption(&value, "jac_recalc_step", mxOptions))
            ekf.system.jac_recalc_step= value;

        if(tryGetOption(&value, "max_steps", mxOptions))
            ekf.system.max_steps= value;
        
        
        double *doflocked= (double*)mxMalloc(EKF::nbrdof*sizeof(double));
        if(tryGetOption(doflocked, "doflocked", mxOptions, EKF::nbrdof, 1)) {
            for(int i= 0; i<EKF::nbrdof; i++)
                ekf.system.doflocked[i]= doflocked[i]!=0.0;
        }
        mxFree(doflocked);
    }
    
    double *q0= mxGetPr(prhs[rhs_idx_q0]);
    double *dq0= mxGetPr(prhs[rhs_idx_dq0]);
    for(int i= 0; i<EKF::nbrdof; ++i) {
        ekf.system.q(i)= q0[i];
        ekf.system.qd(i)= dq0[i];
    }
    
    double *x_ul= mxGetPr(prhs[rhs_idx_x_ul]);
    for(int i=0; i<EKF::nbrstates; ++i)
        ekf.x_ul(i)= x_ul[i];
    double *x_ll= mxGetPr(prhs[rhs_idx_x_ll]);
    for(int i=0; i<EKF::nbrstates; ++i)
        ekf.x_ll(i)= x_ll[i];

    if(nrhs>rhs_idx_Q && !mxIsEmpty(prhs[rhs_idx_Q])) {
        double *Q= mxGetPr(prhs[rhs_idx_Q]);
        for(int i=0; i<EKF::nbrstates; ++i)
            for(int j=0; j<EKF::nbrstates; ++j)
                ekf.ekfQ(i, j)= Q[i + j*EKF::nbrstates];
    }
    if(nrhs>rhs_idx_R && !mxIsEmpty(prhs[rhs_idx_R])) {
        double *R= mxGetPr(prhs[rhs_idx_R]);
        for(int i=0; i<EKF::nbrout; ++i)
            for(int j=0; j<EKF::nbrout; ++j)
                ekf.ekfR(i, j)= R[i + j*EKF::nbrout];
    }
    if(nrhs>rhs_idx_P && !mxIsEmpty(prhs[rhs_idx_P])) {
        double *P= mxGetPr(prhs[rhs_idx_P]);
        for(int i=0; i<EKF::nbrstates; ++i)
            for(int j=0; j<EKF::nbrstates; ++j)
                ekf.ekfSigma(i, j)= P[i + j*EKF::nbrstates];
    }
    if(nrhs>rhs_idx_T_adapt)
        ekf.T_adapt= mxGetScalar(prhs[rhs_idx_T_adapt]);
        
    if(nrhs>rhs_idx_adaptScale && !mxIsEmpty(prhs[rhs_idx_adaptScale])) {
        double *adaptScale= mxGetPr(prhs[rhs_idx_adaptScale]);
        for(int i=0; i<EKF::nbrout; ++i)
            ekf.adaptScale(i)= adaptScale[i];
    }

    if(nrhs>rhs_idx_fixedQxx && !mxIsEmpty(prhs[rhs_idx_fixedQxx])) {
        double *fixedQxx= mxGetPr(prhs[rhs_idx_fixedQxx]);
        for(int i=0; i<EKF::nbrstates; ++i)
            ekf.fixedQxx(i)= fixedQxx[i];
    }

    if(nrhs>rhs_idx_fixedRxx && !mxIsEmpty(prhs[rhs_idx_fixedRxx])) {
        double *fixedRxx= mxGetPr(prhs[rhs_idx_fixedRxx]);
        for(int i=0; i<EKF::nbrout; ++i)
            ekf.fixedRxx(i)= fixedRxx[i];
    }

    double *u= mxGetPr(prhs[rhs_idx_u]);
    double *y_meas= mxGetPr(prhs[rhs_idx_y]);
    
    double ts= mxGetScalar(prhs[rhs_idx_ts]);
    
    ekf.system.precalcConsts();

    plhs[lhs_idx_q]= mxCreateDoubleMatrix(EKF::nbrdof, mxGetN(prhs[rhs_idx_u]), mxREAL);
    plhs[lhs_idx_qd]= mxCreateDoubleMatrix(EKF::nbrdof, mxGetN(prhs[rhs_idx_u]), mxREAL);
    plhs[lhs_idx_qdd]= mxCreateDoubleMatrix(EKF::nbrdof, mxGetN(prhs[rhs_idx_u]), mxREAL);
    plhs[lhs_idx_y]= mxCreateDoubleMatrix(EKF::nbrout, mxGetN(prhs[rhs_idx_u]), mxREAL);

    if(nlhs>lhs_idx_d_norm) {
        plhs[lhs_idx_d_norm]= mxCreateDoubleMatrix(1, mxGetN(prhs[rhs_idx_u]), mxREAL);
    }    
    if(nlhs>lhs_idx_p_xx) {
        plhs[lhs_idx_p_xx]= mxCreateDoubleMatrix(ekf.nbrstates, mxGetN(prhs[rhs_idx_u]), mxREAL);
    }    
    if(nlhs>lhs_idx_r_xx) {
        plhs[lhs_idx_r_xx]= mxCreateDoubleMatrix(ekf.nbrout, mxGetN(prhs[rhs_idx_u]), mxREAL);
    }    
    if(nlhs>lhs_idx_s_xx) {
        plhs[lhs_idx_s_xx]= mxCreateDoubleMatrix(ekf.nbrout, mxGetN(prhs[rhs_idx_u]), mxREAL);
    }    
    
    EKF::VecO y_meas_vec;
    
    std::clock_t startcputime = std::clock();
    try {
        for(int i= 0; i<(int)mxGetN(prhs[rhs_idx_u]); ++i) {
            if(i>0) {
                for(int j= 0; j<ekf.nbrin; ++j)
                    ekf.system.u(j)= u[j + (i-1)*ekf.nbrin];
                
                for(int j= 0; j<ekf.nbrout; ++j)
                    y_meas_vec(j)= y_meas[j + i*ekf.nbrout];
                    
                ekf.next(ts, y_meas_vec);
            }
            for(int j=0; j<ekf.nbrdof; ++j)
                mxGetPr(plhs[lhs_idx_q])[j + i*ekf.nbrdof]= ekf.system.q(j);
            for(int j=0; j<ekf.nbrdof; ++j)
                mxGetPr(plhs[lhs_idx_qd])[j +i*ekf.nbrdof]= ekf.system.qd(j);
            for(int j=0; j<ekf.nbrdof; ++j)
                mxGetPr(plhs[lhs_idx_qdd])[j + i*ekf.nbrdof]= ekf.system.qdd(j);
            for(int j=0; j<ekf.nbrout; ++j)
                mxGetPr(plhs[lhs_idx_y])[j + i*ekf.nbrout]= ekf.system.y(j);
            
            if(nlhs>lhs_idx_d_norm) {
                mxGetPr(plhs[lhs_idx_d_norm])[i]= ekf.d_norm;
            }    
            if(nlhs>lhs_idx_p_xx) {
                for(int j=0; j<ekf.nbrstates; ++j)
                    mxGetPr(plhs[lhs_idx_p_xx])[j + i*ekf.nbrstates]= ekf.ekfSigma_pred(j, j);
            }    
            if(nlhs>lhs_idx_r_xx) {
                for(int j=0; j<ekf.nbrout; ++j)
                    mxGetPr(plhs[lhs_idx_r_xx])[j + i*ekf.nbrout]= ekf.ekfR(j, j);
            }    
            if(nlhs>lhs_idx_s_xx) {
                for(int j=0; j<ekf.nbrout; ++j)
                mxGetPr(plhs[lhs_idx_s_xx])[j + i*ekf.nbrout]= ekf.ekfS(j, j);
            }    
        }
    } catch(std::exception &e) {
        mexWarnMsgIdAndTxt("CADyn:EKF", "Error EKF: %s", e.what());
    }
    double cpu_duration = (std::clock() - startcputime) / (double)CLOCKS_PER_SEC;
    
        
    if(nlhs>lhs_idx_Q) {
        plhs[lhs_idx_Q]= mxCreateDoubleMatrix(EKF::nbrstates, EKF::nbrstates, mxREAL);
        for(int i=0; i<EKF::nbrstates; ++i)
            for(int j=0; j<EKF::nbrstates; ++j)
                mxGetPr(plhs[lhs_idx_Q])[i + j*EKF::nbrstates]= ekf.ekfQ(i, j);
    }    
    if(nlhs>lhs_idx_R) {
        plhs[lhs_idx_R]= mxCreateDoubleMatrix(EKF::nbrout, EKF::nbrout, mxREAL);
        for(int i=0; i<EKF::nbrout; ++i)
            for(int j=0; j<EKF::nbrout; ++j)
                mxGetPr(plhs[lhs_idx_R])[i + j*EKF::nbrout]= ekf.ekfR(i, j);
    }    
    if(nlhs>lhs_idx_P) {
        plhs[lhs_idx_P]= mxCreateDoubleMatrix(EKF::nbrstates, EKF::nbrstates, mxREAL);
        for(int i=0; i<EKF::nbrstates; ++i)
            for(int j=0; j<EKF::nbrstates; ++j)
                mxGetPr(plhs[lhs_idx_P])[i + j*EKF::nbrstates]= ekf.ekfSigma(i, j);
    }    
    if(nlhs>lhs_idx_time) {
        plhs[lhs_idx_time]= mxCreateDoubleMatrix(1, 1, mxREAL);
        mxGetPr(plhs[lhs_idx_time])[0]= cpu_duration;
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
