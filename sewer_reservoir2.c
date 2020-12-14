/*
 Sewer model with five pollutants, flow rate and temperature variation in wastewater and sewer pipe
 *Flow rate is modelled using kinematic wave equation
 *Temperature model considers: i) heat transfer between wastewater and in-sewer air; 
 * ii) convective heat transfer from wastewater to pipe surface; and 
 * iii) conductive heat transfer between wastewater and soil through sewer pipe
 *
 *Ramesh Saagi, IEA, Lund University
 *May 2018
 */

#define S_FUNCTION_NAME  sewer_reservoir2
#define S_FUNCTION_LEVEL 2

#include "simstruc.h"
#include <math.h>

#define XINIT   ssGetSFcnParam(S,0) /* Initial values*/
#define DIM ssGetSFcnParam(S,1)     /* Pipe dimensions and mannings coefficient */
#define PARAM ssGetSFcnParam(S,2)   /* Temperature model parameters*/
#define Aarray   ssGetSFcnParam(S,3) /*A/Afull table to calculate flow rate according to Mannings equation */
#define Siarray   ssGetSFcnParam(S,4)/*Si/Sifull table to calculate flow rate according to Mannings equation */
#define Harray   ssGetSFcnParam(S,5) /*H/Hfull table to determine flow depth */
#define ONOFF   ssGetSFcnParam(S,6) /*ONOFF switch for different heat transfer modes */

static void mdlInitializeSizes(SimStruct *S)
{
    ssSetNumInputPorts(             S,1); 
    ssSetNumOutputPorts(            S,1);   
    ssSetInputPortWidth(            S, 0, 10); /* CODsol, CODpart, NH4-N, TKN, PO4, Ppart,flow rate, water temperature, air temperature and soil temperature*/
    ssSetOutputPortWidth(           S, 0, 19); /* 10 outputs + 8 states*/   
    ssSetNumContStates(             S, 9);   /* number of continuous states - 6 pollutants, Vol, Tww, Tp  */
    ssSetNumDiscStates(             S, 0);  
    ssSetInputPortDirectFeedThrough(S, 0, 1);
    ssSetNumSFcnParams(             S, 7);
    ssSetNumSampleTimes(            S, 1);
    ssSetSimStateCompliance(        S, USE_DEFAULT_SIM_STATE);
    ssSetOptions(                   S,
                                        SS_OPTION_EXCEPTION_FREE_CODE |
                                        SS_OPTION_USE_TLC_WITH_ACCELERATOR);
    
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
    x0[i] = mxGetPr(XINIT)[i];     /* Initial pollutant mass, volume and temperature */
   }
}

#define MDL_UPDATE
static void mdlUpdate(SimStruct *S, int_T tid)
{
 
}
static real_T lookup_h(SimStruct *S)
{
    real_T *x = ssGetContStates(S);
    
    real_T xarray[51],yarray[51]; /*arrays for area and height lookup table */
    const real_T length = mxGetPr(DIM)[0];/*sewer length m*/
    const real_T diameter=mxGetPr(DIM)[1]; /*pipe diameter m*/
    const real_T Afull=0.7854*diameter*diameter; /*Area of fully filled pipe m2 */
    const real_T Hfull=diameter; /*Height of fully filled pipe = pipe diameter m */
    real_T A,H, Aratio, Hratio;
    int i;
    
    for(i=0;i<51;i++)
    {xarray[i]=mxGetPr(Aarray)[i];
    yarray[i]=mxGetPr(Harray)[i];}
            
    /*A=*uPtrs[0];*/
    A=x[6]/length;
    Aratio=A/Afull;        
        
        if (Aratio<=xarray[0])
        {    Hratio=yarray[0];        }
        else if (Aratio >=xarray[50])
        {  Hratio=yarray[50];        }
        else 
        {
            for (i=0;i<50;i++)
            {
              if(Aratio==xarray[i])
              {  Hratio=yarray[i];}
              else if(Aratio>xarray[i] && Aratio<xarray[i+1])
              {Hratio=((Aratio-xarray[i])*(yarray[i+1]-yarray[i])/(xarray[i+1]-xarray[i]))+yarray[i];}
            }
        }
       
    
    H=Hfull*Hratio;
    return H;
}
static real_T lookup_si(SimStruct *S)
{
    real_T *x = ssGetContStates(S);
    
    real_T xarray[51],yarray[51];/*arrays for area and Si value lookup table */
    const real_T length = mxGetPr(DIM)[0];/*sewer length m*/
    const real_T diameter=mxGetPr(DIM)[1];/*pipe diameter m*/
    const real_T Afull=0.7854*diameter*diameter; /*Area of fully filled pipe m2 */
    const real_T Rfull=0.25*diameter; /* Hydraulic radius for fully filled pipe m */
    const real_T Sifull=Afull*pow(Rfull,0.666); /* Si value for fully filled pipe */
    real_T A,Si, Aratio, Siratio;
    int i;
    
    for(i=0;i<51;i++)
    {xarray[i]=mxGetPr(Aarray)[i];
    yarray[i]=mxGetPr(Siarray)[i];}
              
    A=x[6]/length;
    Aratio=A/Afull;        
        
        if (Aratio<=xarray[0])
        {    Siratio=yarray[0];        }
        else if (Aratio >=xarray[50])
        {  Siratio=yarray[50];        }
        else 
        {
            for (i=0;i<50;i++)
            {
              if(Aratio==xarray[i])
              {  Siratio=yarray[i];}
              else if(Aratio>xarray[i] && Aratio<xarray[i+1])
              {Siratio=((Aratio-xarray[i])*(yarray[i+1]-yarray[i])/(xarray[i+1]-xarray[i]))+yarray[i];}
            }
        }
       
    
    Si=Sifull*Siratio;
    return Si;
}

static void mdlOutputs(SimStruct *S, int_T tid)
{
    real_T *y = ssGetOutputPortRealSignal(S,0);
    real_T *x = ssGetContStates(S);
    InputRealPtrsType uPtrs = ssGetInputPortRealSignalPtrs(S,0);
    int i;
    
    real_T factor=86400*pow(mxGetPr(DIM)[2],0.5)/mxGetPr(DIM)[3];
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
    real_T Qout,Wp,Si,Apipe,width, H,qwa_test, qwp_test,qps_p_test,qwp_p_test,Re, Pr, Rh,vel, Dh,alphaw, K,qcod_test,Ktotal_w,Ktotal_p, Tair, Tsoil;
    
    Tair=*uPtrs[8];
    Tsoil=*uPtrs[9];
    Apipe=x[6]/length;
    Si=lookup_si(S);
    H=lookup_h(S);
    Rh=pow(Si/Apipe,1.5);/*Hydraulic radius */
    Wp=Apipe/Rh;
    width = 2*pow(H*(diameter-H),0.5);
    Qout=Si*factor;   
    vel=Qout/(Apipe*86400); /*wastewater velocity */
    Dh=4*Rh;/*Hydraulic diameter */
    Re=vel*Dh/(1.43*pow(10,-6)); /*Reynolds number */
    Pr=(1.43*pow(10,-3))*4181/0.6; /* Prandtl number */
    
    alphaw=0.023*pow(Re,0.8)*pow(Pr,0.33)*0.6/Rh; /*heat transfer coefficient - wastewater to sewer pipe wall */

    Ktotal_w = 1/((qws_onoff*(wt*0.5/kp))+(qwp_conv_onoff/alphaw)); /*total heat transfer coefficient - 1/Ktotal = 1/h1 + 1/h2*/
    Ktotal_p = 1/(qws_onoff*((wt*0.5/kp)+(ds/ks)));
 
    qwa_test=qwa_onoff*(hwa*width*length); /*qwa_onoff*(86400/(1000*4181*x[6]))*(hwa*width*length)*(x[7]-Tair); /* heat flux wastewater - insewer air */
    qwp_test=qwp_conv_onoff*(Ktotal_w)*Wp*length; /*qwp_conv_onoff*(86400/(1000*4181*x[6]))*(Ktotal_w)*Wp*length*(x[7]-x[8]);*/
    qwp_p_test=qwp_conv_onoff*(kp*2/wt)*Wp*length;/*qwp_conv_onoff*(86400/(2400*1000*length*wt*Wp))*(kp*2/wt)*Wp*length*(x[8]-x[7]);*/
    qps_p_test=(Ktotal_p)*Wp*length; /*(86400/(2400*1000*length*wt*Wp))*(Ktotal_p)*Wp*length*(x[8]-Tsoil); */
    qcod_test=qcod_onoff*ecod*rcod*exp(0.078*(x[7]-273.15-20));/*qcod_onoff*(86400/(1000*4181))*ecod*rcod*exp(0.078*(x[7]-273.15-20)); /* heat loss due to COD degradation */
      
    y[0]=K*x[0]; 
    y[1]=K*x[1];
    y[2]=K*x[2];
    y[3]=K*x[3];
    y[4]=K*x[4];
    y[5]=K*x[5];
    y[6]=Qout;
    y[7]=x[7];
    y[8]=Tair;
    y[9]=Tsoil;
    y[10]=x[6]; /*Rh;*/
    y[11]=Rh;
    y[12]=H;
    y[13]=Wp; /*Rhsoil;*/
    y[14]=(qwa_test+qwp_test-qcod_test); /*hwa*width*length;*/
    y[15]=qwp_p_test+qps_p_test; /*Ktotal_w;*/
    y[16]=qwp_p_test; /*Ktotal_p;  */
    y[17]=qps_p_test;/*H;/*x[6]; /*x[8];*/
    y[18]=x[8]; /*qcod_test; /*theta;*/
   }

#define MDL_DERIVATIVES 
static void mdlDerivatives(SimStruct *S)
  {
    real_T            *dx   = ssGetdX(S);    
    real_T            *x    = ssGetContStates(S);
    InputRealPtrsType uPtrs = ssGetInputPortRealSignalPtrs(S,0);
    int i;
    
    real_T factor=86400*pow(mxGetPr(DIM)[2],0.5)/mxGetPr(DIM)[3];
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
    real_T Qout,Wp,Si,Apipe,width, H,qwa, qwp,qps_p,qwp_p,Re, Pr, Rh,vel, Dh,alphaw, K,qcod,Ktotal_w,Ktotal_p, Tair, Tsoil;
    real_T theta,Rpipemid,Apipemid,Wpipemid,Rhpipemid,Rpipeend,Apipeend,Wpipeend, Rhpipeend,Rsoil,Wsoil,Asoil,Rhsoil,Rtotal_w,Rtotal_p;
    
    Tair=*uPtrs[8];
    Tsoil=*uPtrs[9];
    Apipe=x[6]/length;
    Si=lookup_si(S);
    H=lookup_h(S);
    Rh=pow(Si/Apipe,1.5);/*Hydraulic radius */
    Wp=Apipe/Rh;
    width = 2*pow(H*(diameter-H),0.5);
    Qout=Si*factor;   
    vel=Qout/(Apipe*86400); /*wastewater velocity */
    Dh=4*Rh;/*Hydraulic diameter */
    Re=vel*Dh/(1.43*pow(10,-6)); /*Reynolds number */
    Pr=(1.43*pow(10,-3))*4181/0.6; /* Prandtl number */
    
   alphaw=0.023*pow(Re,0.8)*pow(Pr,0.33)*0.6/Rh; /*heat transfer coefficient - wastewater to sewer pipe wall */
   Ktotal_w = 1/((qws_onoff*(wt*0.5/kp))+(qwp_conv_onoff/alphaw)); /*total heat transfer coefficient - 1/Ktotal = 1/h1 + 1/h2 + 1/h3*/
   Ktotal_p = 1/(qws_onoff*((wt*0.5/kp)+(ds/ks)));
 
   qwa=qwa_onoff*(86400/(1000*4181*x[6]))*(hwa*width*length*(x[7]-Tair)); /* heat flux wastewater - insewer air */
   qwp=qwp_conv_onoff*(86400/(1000*4181*x[6]))*(Ktotal_w)*Wp*length*(x[7]-x[8]);
   qwp_p=qwp_p=qwp_conv_onoff*(86400/(2400*1000*length*wt*Wp))*(kp*2/wt)*Wp*length*(x[8]-x[7]);
   qps_p=(86400/(2400*1000*length*wt*Wp))*(Ktotal_p)*Wp*length*(x[8]-Tsoil); 
   qcod=qcod_onoff*(86400/(1000*4181))*ecod*rcod*exp(0.078*(x[7]-273.15-20)); /* heat loss due to COD degradation */
    
    K=Qout/x[6]; /*1/Residence time (/d ) */          
    
    dx[0]=*uPtrs[0]-K*x[0]-rcod*x[6];
    for (i=1;i<6; i++)
      {
        dx[i]=*uPtrs[i]-K*x[i];
      }
      dx[6]=*uPtrs[6]-Qout; /*volume balance */
      dx[7]= (*uPtrs[6]*(*uPtrs[7]-x[7])/x[6])-qwa-qwp+qcod; /*energy balance */
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
