%% Load rawdata for Linkoping. 
load('processeddata_malmo.mat')
%T_airtemp1=T_airtemp([1:18989,19000:end],:); % this is to remove extra hour due to change from day light savings time 

%% Process raw data - convert to timetable format for easy postprocessing and also convert units to m3/d and K

%Convert units - flow rate (m3/h to m3/d)and temperature (C to K)
T_down_5min.T=T_down_5min.T+273.15;
T_air_5min.T=T_air_5min.T+273.15;
T_up_5min.T=T_up_5min.T+273.15;

%select the time interval where all the data is available - March 08th to June 29th
startdate=datetime(2019,03,09,00,00,00);
enddate=datetime(2019,06,28,00,00,00);

ind1=find(Q_down_5min.date>=startdate,1,'first');
ind2=find(Q_down_5min.date<=enddate,1,'last');
Q_down_5min=Q_down_5min(ind1:ind2,:);

ind1=find(T_down_5min.date>=startdate,1,'first');
ind2=find(T_down_5min.date<=enddate,1,'last');
T_down_5min=T_down_5min(ind1:ind2,:);

ind1=find(T_air_5min.date>=startdate,1,'first');
ind2=find(T_air_5min.date<=enddate,1,'last');
T_air_5min=T_air_5min(ind1:ind2,:);

ind1=find(T_up_5min.date>=startdate,1,'first');
ind2=find(T_up_5min.date<=enddate,1,'last');
T_up_5min=T_up_5min(ind1:ind2,:);

in_Malmo_up.time=[];
in_Malmo_up.signals.values(:,7)=Q_down_5min.Q;%Downstream flow rate used as there is no upstream data
in_Malmo_up.signals.values(:,8)=T_up_5min.T;
in_Malmo_up.signals.values(:,1:6)=0;

in_airdata.time=[];
in_airdata.signals.values=T_air_5min.T;

td=Q_down_5min.date; %Time variable in date format for plotting
tt=(0:5/(24*60):days(td(end)-td(1)))'; %Time variable in days for plotting