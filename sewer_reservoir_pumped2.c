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

#define S_FUNCTION_NAME  sewer_reservoir_pumped2
#define S_FUNCTION_LEVEL 2

#include "simstruc.h"
#include <math.h>

#define XINIT   ssGetSFcnParam(S,0) /* Initial values*/
#define DIM ssGetSFcnParam(S,1)     /* Pipe dimensions and mannings coefficient */
#define PARAM ssGetSFcnParam(S,2)   /* Temperature model parameters*/
#define ONOFF   ssGetSFcnParam(S,3) /*ONOFF switch for different heat transfer modes */

static void mdlInitializeSizes(SimStruct *S)
{
    ssSetNumInputPorts(             S,1); 
    ssSetNumOutputPorts(            S,1);   
    ssSetInputPortWidth(            S, 0, 10); /* CODsol, CODpart, NH4-N, TKN, PO4, Ppart,flow rate, water temperature, air temp and soil temp*/
    ssSetOutputPortWidth(           S, 0, 14); /* 14 outputs */   
    ssSetNumContStates(             S, 9);   /* number of continuous states - 6 pollutants, vol, T  */
    ssSetNumDiscStates(             S, 0);  
    ssSetInputPortDirectFeedThrough(S, 0, 1);
    ssSetNumSFcnParams(             S, 4);
    ssSetNumSampleTimes(            S, 1);
    ssSetSimStateCompliance(        S, USE_DEFAULT_SIM_STATE);
    ssSetOptions(                   S,
                                        SS_OPTION_EXCEPTION_FREE_CODE |
                                        SS_OPTION_USE_TLC_WITH_ACCELERATOR); /*Removed SS_OPTION_WORKS_WITH_CODE_REUSE | to avoid warnings */
    
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
    for (i=0;i<9; i++)
   {
    x0[i] = mxGetPr(XINIT)[i];     /* Initial pollutant conc, volume and temperature */
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
    
    const real_T length = mxGetPr(DIM)[0];/*m*/
    const real_T diameter=mxGetPr(DIM)[1]; /*m*/
    const real_T hwa=mxGetPr(PARAM)[0];
    const real_T kp=mxGetPr(PARAM)[1];
    const real_T ks=mxGetPr(PARAM)[2];
    const real_T wt=mxGetPr(PARAM)[3];
    const real_T ds=mxGetPr(PARAM)[4];
    const real_T ecod=mxGetPr(PARAM)[5];/*   */
    const real_T rcod=mxGetPr(PARAM)[6];/**/
    const real_T qwa_onoff=mxGetPr(ONOFF)[0];/**/
    const real_T qws_onoff=mxGetPr(ONOFF)[1];/**/
    const real_T qwp_conv_onoff=mxGetPr(ONOFF)[2];/**/
    const real_T qcod_onoff=mxGetPr(ONOFF)[3];/**/
    real_T Qout,Wp,Apipe,H,qwa, qwp_test,qwp_p_test,qps_p_test,Re, Pr, Rh,vel, Dh,alphaw, qcod_test, Ktotal_w,Ktotal_p,vol,Tair, Tsoil;
    
    Tair=*uPtrs[8];
    Tsoil=*uPtrs[9];
    Apipe=0.7854*diameter*diameter; /*Area of fully filled pipe m2 */
    Rh=0.25*diameter; /* Hydraulic radius for fully filled pipe m */
    vol=Apipe*length;
    H=diameter;
    Wp=Apipe/Rh;
    Dh=4*Rh;/*Hydraulic diameter */
    Qout=*uPtrs[6]; /*pumped flow rate */
    vel=Qout/(Apipe*86400); /*wastewater velocity */
    Re=vel*Dh/(1.43*pow(10,-6)); /*Reynolds number */
    Pr=(1.43*pow(10,-3))*4181/0.6; /* Prandtl number */
 
    
    alphaw=0.023*pow(Re,0.8)*pow(Pr,0.33)*0.6/Rh; /*heat transfer coefficient - wastewater to sewer pipe wall */
    Ktotal_w =  1/((qws_onoff*(wt*0.5/kp))+(qwp_conv_onoff/alphaw)); /*total heat transfer coefficient - 1/Ktotal = 1/h1 + 1/h2 + 1/h3*/
    Ktotal_p = 1/(qws_onoff*((wt*0.5/kp)+(ds/ks)));

    qwa=0; /*qwa_onoff*(86400/(1000*4181*x[6]))*(hwa*width*length*(x[7]-Tair)); /* heat flux wastewater - insewer air */
    qwp_test=qwp_conv_onoff*(Ktotal_w)*Wp*length;/*qwp_conv_onoff*(86400/(1000*4181*x[6]))*(Ktotal_w)*Wp*length*(x[7]-x[8]);*/
    qwp_p_test=qwp_conv_onoff*(kp*2/wt)*Wp*length;
    qps_p_test=(Ktotal_p)*Wp*length;
    qcod_test=qcod_onoff*ecod*rcod*exp(0.078*(x[7]-273.15-20)); /* heat loss due to COD degradation */
    
    
         
    y[0]=x[0]; 
    y[1]=x[1];
    y[2]=x[2];
    y[3]=x[3];
    y[4]=x[4];
    y[5]=x[5];
    y[6]=Qout;
    y[7]=x[7];
    y[8]=Tair;
    y[9]=Tsoil;
    y[10]=qwp_test-qcod_test; /*x[8];*/
    y[11]=qwp_p_test+qps_p_test; /*Ktotal_w;*/
    y[12]=qps_p_test; /*Ktotal_p;*/
    y[13]=qcod_test; /*Rsoil;*/

}

#define MDL_DERIVATIVES 
static void mdlDerivatives(SimStruct *S)
  {
    real_T            *dx   = ssGetdX(S);    
    real_T            *x    = ssGetContStates(S);
    InputRealPtrsType uPtrs = ssGetInputPortRealSignalPtrs(S,0);
    int i;
    
    const real_T length = mxGetPr(DIM)[0];/*m*/
    const real_T diameter=mxGetPr(DIM)[1]; /*m*/
    const real_T hwa=mxGetPr(PARAM)[0];
    const real_T kp=mxGetPr(PARAM)[1];
    const real_T ks=mxGetPr(PARAM)[2];
    const real_T wt=mxGetPr(PARAM)[3];
    const real_T ds=mxGetPr(PARAM)[4];
    const real_T ecod=mxGetPr(PARAM)[5];/*   */
    const real_T rcod=mxGetPr(PARAM)[6];/**/
    const real_T qwa_onoff=mxGetPr(ONOFF)[0];/**/
    const real_T qws_onoff=mxGetPr(ONOFF)[1];/**/
    const real_T qwp_conv_onoff=mxGetPr(ONOFF)[2];/**/
    const real_T qcod_onoff=mxGetPr(ONOFF)[3];/**/
    real_T Qout,Wp,Apipe,H,qwa, qwp,qwp_p,qps_p,Re, Pr, Rh,vel, Dh,alphaw, qcod, Ktotal_w,Ktotal_p,vol,Tair, Tsoil;
    
    Tair=*uPtrs[8];
    Tsoil=*uPtrs[9];
    Apipe=0.7854*diameter*diameter; /*Area of fully filled pipe m2 */
    Rh=0.25*diameter; /* Hydraulic radius for fully filled pipe m */
    vol=Apipe*length;
    H=diameter;
    Wp=Apipe/Rh;
    Dh=4*Rh;/*Hydraulic diameter */
    Qout=*uPtrs[6]; /*pumped flow rate */
    vel=Qout/(Apipe*86400); /*wastewater velocity */
    Re=vel*Dh/(1.43*pow(10,-6)); /*Reynolds number */
    Pr=(1.43*pow(10,-3))*4181/0.6; /* Prandtl number */
    
    alphaw=0.023*pow(Re,0.8)*pow(Pr,0.33)*0.6/Rh; /*heat transfer coefficient - wastewater to sewer pipe wall */
    Ktotal_w =  1/((qws_onoff*(wt*0.5/kp))+(qwp_conv_onoff/alphaw)); /*total heat transfer coefficient - 1/Ktotal = 1/h1 + 1/h2 + 1/h3*/
    Ktotal_p = 1/(qws_onoff*((wt*0.5/kp)+(ds/ks)));

    qwa=0; /*qwa_onoff*(86400/(1000*4181*x[6]))*(hwa*width*length*(x[7]-Tair)); /* heat flux wastewater - insewer air */
    qwp=qwp_conv_onoff*(86400/(1000*4181*x[6]))*(Ktotal_w)*Wp*length*(x[7]-x[8]);
    qwp_p=qwp_conv_onoff*(86400/(2400*1000*length*wt*Wp))*(kp*2/wt)*Wp*length*(x[8]-x[7]);
    qps_p=(86400/(2400*1000*length*wt*Wp))*(Ktotal_p)*Wp*length*(x[8]-Tsoil);
    qcod=qcod_onoff*(86400/(1000*4181))*ecod*rcod*exp(0.078*(x[7]-273.15-20)); /* heat loss due to COD degradation */
    
    dx[0]=((Qout/vol)*(*uPtrs[0]-x[0]))-rcod*x[6];
    for (i=1;i<6; i++)
      {
        dx[i]=(Qout/vol)*(*uPtrs[i]-x[i]);
      }
    dx[6]=0; /*volume balance */
    dx[7]= (Qout/vol)*(*uPtrs[7]-x[7])-qwa-qwp+qcod; /*energy balance */
    dx[8]= -qwp_p-qps_p; /*energy balance */

  }

static void mdlTerminate(SimStruct *S)
{
}

#ifdef  MATLAB_MEX_FILE    /* Is this file being compiled as a MEX-file? */
#include "simulink.c"      /* MEX-file interface mechanism */
#else
#include "cg_sfun.h"       /* Code generation registration function */
#endif
