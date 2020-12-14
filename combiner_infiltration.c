/*
 * Flow combiner for infiltration to sewers
 * Ramesh Saagi, IEA, Lund University */

#define S_FUNCTION_NAME combiner_infiltration

#include "simstruc.h"
#include <math.h>

 /*
 * mdlInitializeSizes - initialize the sizes array
 */
static void mdlInitializeSizes(SimStruct *S)
{
    ssSetNumContStates(    S, 0);   /* number of continuous states                  (CODs,CODp, SNH, TKN,PO4, Ppart,Q,T)       */
    ssSetNumDiscStates(    S, 0);   /* number of discrete states                                                        */
    ssSetNumInputs(        S, 20);  /* number of inputs                              (CODs,CODp, SNH, TKN,PO4, Ppart,Q,T)      */
    ssSetNumOutputs(       S, 8);   /* number of outputs                                                                */
    ssSetDirectFeedThrough(S, 1);   /* direct feedthrough flag                                                          */
    ssSetNumSampleTimes(   S, 1);   /* number of sample times                                                           */
    ssSetNumSFcnParams(    S, 0);   /* number of input arguments                                                        */
    ssSetNumRWork(         S, 0);   /* number of real work vector elements                                              */
    ssSetNumIWork(         S, 0);   /* number of integer work vector elements                                           */
    ssSetNumPWork(         S, 0);   /* number of pointer work vector elements                                           */
}

/*
 * mdlInitializeSampleTimes - initialize the sample times array
 */
static void mdlInitializeSampleTimes(SimStruct *S)
{
    ssSetSampleTime(S, 0, CONTINUOUS_SAMPLE_TIME);
    ssSetOffsetTime(S, 0, 0.0);
}
/*
 * mdlInitializeConditions - initialize the states
 */
static void mdlInitializeConditions(double *x0, SimStruct *S)
{
}
/*
 * mdlOutputs - compute the outputs
 */
static void mdlOutputs(double *y, double *x, double *u, SimStruct *S, int tid)
{
int i;

for (i=0;i<7; i++)
   {
    y[i] = (u[i]+u[8+i]);     /* Output pollutant load and flow rate */
   }
if (u[6]+u[14]> 0.1)
{y[7]=(u[7]*u[6]+u[15]*u[14])/(u[6]+u[14]);}
else
{y[7]=20+273.15;}           /* Special case when input flow rate and infiltration are both zero, 20C is assumed as the temp now. In future, ideal to use temp from the preceeding input data. */ 

}  
/*
 * mdlUpdate - perform action at major integration time step
 */
static void mdlUpdate(double *x, double *u, SimStruct *S, int tid)
{
}
/*
 * mdlDerivatives - compute the derivatives
 */
static void mdlDerivatives(double *dx, double *x, double *u, SimStruct *S, int tid)
{
}
/*
 * mdlTerminate - called when the simulation is terminated.
 */
static void mdlTerminate(SimStruct *S)
{
}

#ifdef	MATLAB_MEX_FILE    /* Is this file being compiled as a MEX-file? */
#include "simulink.c"      /* MEX-file interface mechanism */
#else
#include "cg_sfun.h"       /* Code generation registration function */
#endif
