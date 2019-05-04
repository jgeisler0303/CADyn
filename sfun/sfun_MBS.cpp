#define makeMBSystemSFunName(c) sfun_ ## c
#define exp_makeMBSystemSFunName(c) makeMBSystemSFunName(c)

#define S_FUNCTION_NAME exp_makeMBSystemSFunName(SFUNMBSystem)
#define S_FUNCTION_LEVEL 2

#include "simstruc.h"
#include <math.h> 

#define makeMBSystemClass(c) c ## System
#define exp_makeMBSystemClass(c) makeMBSystemClass(c)
#define MBSystemClass exp_makeMBSystemClass(SFUNMBSystem)

#define stringify(x)  #x
#define expand_and_stringify(x) stringify(x ## System.hpp)
#define INCL_FILE_STR(x) expand_and_stringify(x)
#include INCL_FILE_STR(SFUNMBSystem)

#define MAX_ERR_MSG_LEN 256


#define NPARAMS 3

#define PARAM_TS_IDX 0
#define PARAM_MX_TS ssGetSFcnParam(S, PARAM_TS_IDX)
#define PARAM_P_TS mxGetPr(PARAM_MX_TS)
#define PARAM_TS mxGetScalar(PARAM_MX_TS)

#define PARAM_IC_Q_IDX 1
#define PARAM_MX_IC_Q ssGetSFcnParam(S, PARAM_IC_Q_IDX)
#define PARAM_P_IC_Q mxGetPr(PARAM_MX_IC_Q)

#define PARAM_IC_QD_IDX 2
#define PARAM_MX_IC_QD ssGetSFcnParam(S, PARAM_IC_QD_IDX)
#define PARAM_P_IC_QD mxGetPr(PARAM_MX_IC_QD)



#ifdef MATLAB_MEX_FILE
#define MDL_CHECK_PARAMETERS
static void mdlCheckParameters(SimStruct *S) {
    static char msg[MAX_ERR_MSG_LEN];
    
    if(!(mxIsDouble(PARAM_MX_TS) && mxIsScalar(PARAM_MX_TS))) {
         ssSetErrorStatus(S, "Parameter no 1 ('ts') must be a real scalar");
         return;
    }
    if(PARAM_TS<0.0 || !mxIsFinite(PARAM_TS)) {
         ssSetErrorStatus(S, "Parameter no 1 ('ts') must be a positive");
         return;
    }
    
    if(!(mxIsDouble(PARAM_MX_IC_Q) && !mxIsSparse(PARAM_MX_IC_Q) && mxGetNumberOfElements(PARAM_MX_IC_Q)==MBSystemClass::nStates)) {
        snprintf(msg, MAX_ERR_MSG_LEN-1, "Parameter no 2 ('ic_q') must be a real dense vector with length %d", MBSystemClass::nStates);
        ssSetErrorStatus(S, msg);
        return;
    }
    
    if(!(mxIsDouble(PARAM_MX_IC_QD) && !mxIsSparse(PARAM_MX_IC_QD) && mxGetNumberOfElements(PARAM_MX_IC_QD)==MBSystemClass::nStates)) {
        snprintf(msg, MAX_ERR_MSG_LEN-1, "Parameter no 3 ('ic_qd') must be a real dense vector with length %d", MBSystemClass::nStates);
        ssSetErrorStatus(S, msg);
        return;
    }
}
#endif

static void mdlInitializeSizes(SimStruct *S) {
    ssSetNumSFcnParams(S, NPARAMS);
#if defined(MATLAB_MEX_FILE)
    if (ssGetNumSFcnParams(S) != ssGetSFcnParamsCount(S)) return;
    mdlCheckParameters(S);
    if (ssGetErrorStatus(S) != NULL) return;
#endif
    
    for(int i= 0; i<NPARAMS; i++) ssSetSFcnParamTunable(S, i, false);

    ssSetNumContStates(S, 0);
    ssSetNumDiscStates(S, 0);

    ssSetNumInputPorts(S, MBSystemClass::nInputs);
    for(int i= 0; i<MBSystemClass::nInputs; i++) {
        ssSetInputPortWidth(S, i, 1);
        ssSetInputPortOptimOpts(S, i, SS_NOT_REUSABLE_AND_GLOBAL);
        ssSetInputPortDirectFeedThrough(S, i, true);
    }
    

    ssSetNumOutputPorts(S, 3);
    for(int i= 0; i<3; i++) {
        ssSetOutputPortWidth(S, i, MBSystemClass::nStates);
        ssSetOutputPortOptimOpts(S, i, SS_NOT_REUSABLE_AND_GLOBAL);
    }
    

    ssSetNumSampleTimes(S, 1);
    ssSetNumRWork(S, 1);
    ssSetNumIWork(S, 0);
    ssSetNumPWork(S, 1);
    ssSetNumModes(S, 0);
    ssSetNumNonsampledZCs(S, 0);

    ssSetOptions(S, SS_OPTION_EXCEPTION_FREE_CODE);
}

static void mdlInitializeSampleTimes(SimStruct *S) {
    ssSetSampleTime(S, 0, PARAM_TS);
    ssSetOffsetTime(S, 0, 0.0);
}

#define MDL_START
void mdlStart(SimStruct *S) {
    MBSystemClass *sys= new MBSystemClass();
    ssGetPWorkValue(S, 0)= sys;
    
    for(int i= 0; i<MBSystemClass::nStates; i++) sys->q(i)= PARAM_P_IC_Q[i];
    for(int i= 0; i<MBSystemClass::nStates; i++) sys->qd(i)= PARAM_P_IC_QD[i];
    sys->qdd.setZero();
    
    ssSetRWorkValue(S, 0, PARAM_TS);
}

static void mdlOutputs(SimStruct *S, int_T tid) {
    MBSystemClass *sys= static_cast<MBSystemClass *>(ssGetPWorkValue(S, 0));
    
    for(int i= 0; i<MBSystemClass::nInputs; i++)
        sys->u[i]= (*ssGetInputPortRealSignalPtrs(S, i)[0]);

    if(ssIsFirstInitCond(S)) {
        sys->t= ssGetT(S);
        double dummy;
        if(sys->newmarkOneStep(0.0, dummy)!=0) {
            ssSetErrorStatus(S, "Error initializing Newmark integrator for system ");
            return;
        }
    } else {
        if(!sys->newmarkInterval(ssGetT(S), *ssGetRWork(S), PARAM_TS)) {
            ssSetErrorStatus(S, "Error in Newmark integrator for system ");
            return;
        }
    }
    
    for(int i= 0; i<MBSystemClass::nStates; i++) ssGetOutputPortRealSignal(S, 0)[i]= sys->q[i];
    for(int i= 0; i<MBSystemClass::nStates; i++) ssGetOutputPortRealSignal(S, 1)[i]= sys->qd[i];
    for(int i= 0; i<MBSystemClass::nStates; i++) ssGetOutputPortRealSignal(S, 2)[i]= sys->qdd[i];
}

static void mdlTerminate(SimStruct *S) {
    MBSystemClass *sys= static_cast<MBSystemClass *>(ssGetPWorkValue(S, 0));
    
    delete sys;
    sys= nullptr;
}

#ifdef  MATLAB_MEX_FILE    /* Is this file being compiled as a MEX-file? */
#include "simulink.c"      /* MEX-file interface mechanism */
#else
#include "cg_sfun.h"       /* Code generation registration function */
#endif

