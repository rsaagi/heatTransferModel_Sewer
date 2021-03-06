/*
 Sewer model with five pollutants, flow rate and temperature variation
 *Flow rate is modelled using kinematic wave equation
 *Temperature model considers: i) heat transfer between wastewater and in-sewer air; 
 * ii) convective heat transfer from wastewater to pipe surface; and 
 * iii) conductive heat transfer between wastewater and soil through sewer pipe
 *
 *Ramesh Saagi, IEA, Lund University
 *May 2018
 */

#define S_FUNCTION_NAME  sewer_reservoir_lumped
#define S_FUNCTION_LEVEL 2

#include "simstruc.h"
#include <math.h>

#define XINIT   ssGetSFcnParam(S,0) /* Initial values*/
#define PARAM ssGetSFcnParam(S,1)   /* Temperature model parameters*/
#define QAVG ssGetSFcnParam(S,2)   /* Temperature model parameters*/

static void mdlInitializeSizes(SimStruct *S)
{
    ssSetNumInputPorts(             S,1); 
    ssSetNumOutputPorts(            S,1);   
    ssSetInputPortWidth(            S, 0, 9); /* CODsol, CODpart, NH4-N, TKN, PO4, Ppart,flow rate, water temperature, air temperature*/
    ssSetOutputPortWidth(           S, 0, 9); /* 9 outputs*/   
    ssSetNumContStates(             S, 8);   /* number of continuous states - 6 pollutants, Vol, T  */
    ssSetNumDiscStates(             S, 0);  
    ssSetInputPortDirectFeedThrough(S, 0, 1);
    ssSetNumSFcnParams(             S, 3);
    ssSetNumSampleTimes(            S, 1);
    ssSetSimStateCompliance(        S, USE_DEFAULT_SIM_STATE);
    ssSetOptions(                   S,
                                        SS_OPTION_EXCEPTION_FREE_CODE |
                                        SS_OPTION_USE_TLC_WITH_ACCELERATOR); /*Removed SS_OPTION_WORKS_WITH_CODE_REUSE | as it throws a warning*/
    
    if (ssGetNumSFcnParams(S) != ssGetSFcnParamsCount(S)) return; /* Parameter mismatch will be reported by Simulink */
}

static void mdlInitializeSampleTimes(SimStruct *S)
{
    ssSetSampleTime(S, 0, CONTINUOUS_SAMPLE_TIME);
    ssSetOffsetTime(S, 0, 0.0);
    ssSetModelReferenceSampleTimeDefaultInheritance(S); 
}

#define MDL_INITIALIZE_CONDITIONS
static void mdlInitializeConditions(SimStruct *S)
{
real_T *x0 = ssGetContStates(S);
int i;
    for (i=0;i<8; i++)
   {
    x0[i] = mxGetPr(XINIT)[i];     /* Initial pollutant mass, volume and temperature */
   }
}

#define MDL_UPDATE
static void mdlUpdate(SimStruct *S, int_T tid)
{
 
}

static void mdlOutputs(SimStruct *S, int_T tid)
{
    real_T *y = ssGetOutputPortRealSignal(S,0);
    real_T *x = ssGetContStates(S);
    InputRealPtrsType uPtrs = ssGetInputPortRealSignalPtrs(S,0);
      
    const real_T K=mxGetPr(PARAM)[2];
    real_T Tair;
    
    Tair = *uPtrs[8];
   

    y[0]=K*x[0]; 
    y[1]=K*x[1];
    y[2]=K*x[2];
    y[3]=K*x[3];
    y[4]=K*x[4];
    y[5]=K*x[5];
    y[6]=K*x[6]; /*Qout;*/
    y[7]=x[7];
    y[8]=Tair;
   }

#define MDL_DERIVATIVES 
static void mdlDerivatives(SimStruct *S)
  {
    real_T            *dx   = ssGetdX(S);    
    real_T            *x    = ssGetContStates(S);
    InputRealPtrsType uPtrs = ssGetInputPortRealSignalPtrs(S,0);
    int i;
    
    const real_T hmax=mxGetPr(PARAM)[0];
    const real_T n_kflow=mxGetPr(PARAM)[1];
    const real_T K=mxGetPr(PARAM)[2];
    const real_T Qavg=mxGetPr(QAVG)[0];
    
    real_T Qout, Tair, hcalc, kflow;
   
    kflow=Qavg/n_kflow;
    Tair=*uPtrs[8];
    Qout=K*x[6];
    hcalc=((hmax*Qout)/(kflow+Qout));
    
    for (i=0;i<6; i++)
      {
        dx[i]=*uPtrs[i]-K*x[i];
      }
    dx[6]=*uPtrs[6]-Qout; /*volume balance */
    dx[7]= (*uPtrs[6]*(*uPtrs[7]-x[7])/x[6])-((86400/(1000*4181*x[6]))*hcalc*(x[7]-Tair));/*qwa-qws-qcod; /*energy balance */
  }

static void mdlTerminate(SimStruct *S)
{
}

#ifdef  MATLAB_MEX_FILE    /* Is this file being compiled as a MEX-file? */
#include "simulink.c"      /* MEX-file interface mechanism */
#else
#include "cg_sfun.h"       /* Code generation registration function */
#endif
