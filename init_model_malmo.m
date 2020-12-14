% Initialization for sewer model with pollutants, flow rate and temperature
% Ramesh Saagi, IEA, Lund University
% May 2018

%% Input mean values

% Normalized profiles
processrawdata_malmosewer
%Other temperature inputs
Tair=11+273.15; % K
Tsoil=8+273.15; %13.1+273.15; % K
%% Pipe dimensions

%% gravity sewer pipe - p2
length=340; % pipe length for each section (m)
diameter=0.4; % pipe diameter (m)
slope=0.006; % horizontal slope
mannings=0.01; % mannings coefficient for concrete pipe
DIM2=[length, diameter, slope, mannings]; 
%% gravity sewer pipe - p3
length=50; % pipe length for each section (m)
diameter=0.4; % pipe diameter (m)
slope=0.006; % horizontal slope
mannings=0.01; % mannings coefficient for concrete pipe
DIM3=[length, diameter, slope, mannings]; 
%% gravity sewer pipe - p4
length=115; % pipe length for each section (m)
diameter=0.4; % pipe diameter (m)
slope=0.003; % horizontal slope
mannings=0.01; % mannings coefficient for concrete pipe
DIM4=[length, diameter, slope, mannings]; 
%% gravity sewer pipe - p5
length=95; % pipe length for each section (m)
diameter=0.4; % pipe diameter (m)
slope=0.004; % horizontal slope
mannings=0.01; % mannings coefficient for concrete pipe
DIM5=[length, diameter, slope, mannings]; 
%% gravity sewer pipe - p6
length=260; % pipe length for each section (m)
diameter=0.4; % pipe diameter (m)
slope=0.006; % horizontal slope
mannings=0.01; % mannings coefficient for concrete pipe
DIM6=[length, diameter, slope, mannings]; 
%% gravity sewer pipe - p7
length=130; % pipe length for each section (m)
diameter=0.4; % pipe diameter (m)
slope=0.004; % horizontal slope
mannings=0.01; % mannings coefficient for concrete pipe
DIM7=[length, diameter, slope, mannings]; 
%% gravity sewer pipe - p8
length=100; % pipe length for each section (m)
diameter=0.6; % pipe diameter (m)
slope=0.003; % horizontal slope
mannings=0.01; % mannings coefficient for concrete pipe
DIM8=[length, diameter, slope, mannings]; 
%% gravity sewer pipe - p9
length=140; % pipe length for each section (m)
diameter=0.4; % pipe diameter (m)
slope=0.003; % horizontal slope
mannings=0.01; % mannings coefficient for concrete pipe
DIM9=[length, diameter, slope, mannings]; 
%% gravity sewer pipe - p10
length=165; % pipe length for each section (m)
diameter=0.4; % pipe diameter (m)
slope=0.002; % horizontal slope
mannings=0.01; % mannings coefficient for concrete pipe
DIM10=[length, diameter, slope, mannings]; 
%% gravity sewer pipe - p11
length=265; % pipe length for each section (m)
diameter=0.4; % pipe diameter (m)
slope=0.004; % horizontal slope
mannings=0.01; % mannings coefficient for concrete pipe
DIM11=[length, diameter, slope, mannings];

%% Infiltration to sewers
% inf_onoff = 1; %Turn infiltration off. This infiltration is spread equally to all the reservoirs
% noise_onoff=1;
% gwbias=Qin/3; %4850; % Average yearly infiltration 
% amp=0.25; % amplitude fraction for each sub catchment
% freq=2*pi/364;
% phase=-pi*15/24;
% noise_std=0.2; %Noise standard deviation as a fraction of total load 
% gwtemp=Tsoil;
% nreservoirs=10; %No. of reservoirs. This is to spread the infiltration flow equaly to all the reservoirs


%% Heat transfer parameters
load('lookup_ratios.mat') % Lookup tables for Area, Height and Si values based on SWMM hydraulics manual

hwa=5;% 15;%25; %5;%heat transfer coefficient wastewater to in-sewer air W/m2.k
kp=1.3; %5.3; %0.5; %2.3;% %Thermal condictivity of pipe W/m.k
ks=1.5; %5.5; %1.5;%5.5; %5.5;%2.5;% Thermal condictivity of soil W/m.k
wt=0.06; %.04; %0.14; %Pipe thickness m
ds=0.2; %0.1; %soil depth for heat transfer m
ecod=14e6;% reaction enthalpy for COD degradation J/kgCOD
rcod=1e-6; %COD degradation rate in sewers (Huisman et al., 2004) kgCOD/m3.s

qwa_onoff=1; %wastewater - in-sewer air. 1 - ON; 0 - OFF
qws_onoff=1; % wastewater -soil
qwp_conv_onoff=1;%forced convection wastewater due to flow rate
qcod_onoff=1;%cod degradation

%% Create parameter vectors
%input_mean=[CODsol CODpart NH4 TKN PO4 Ppart Qin Tin];

% if(exist('stateset_init.mat','file'))
%     load stateset_init.mat
%     XINIT=stateset_init;
% else
XINIT=[2, 10, 0.5, 0, 0.15, 0.15, 22.85, 7+273.15, 7+273.15]; %initial mass of pollutants (kg), volume (m3) and temperature (K) in sewer pipes
%XINIT_pumped=[2, 10, 0.5, 0, 0.15, 0.15, pi*DIM1(2)*DIM1(2)*DIM1(1)/4, 7+273.15, 7+273.15]; %initial mass of pollutants (kg), volume (m3) and temperature (K) in sewer pipes
PARAM=[hwa, kp, ks, wt, ds, Tair, Tsoil, ecod, rcod]; 
ONOFF=[qwa_onoff, qws_onoff, qwp_conv_onoff, qcod_onoff];

hmax=35000; %50; %1500;%20;
Kres=200; %500; %200;
Qavg=3300; %m3/d
n_flow=0.2;

%PARAM_lumped=[hmax, kflow, Kres]; 
PARAM_lumped=[hmax, n_flow, Kres]; 

clear length

