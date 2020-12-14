% Run sewer mdoel
% Ramesh Saagi, IEA, Lund University
% Oct 2018
%% Initialize
init_model_malmo
tend=days(enddate-startdate);

%% Run the model

tic;sim('sewermodel_malmo_conceptual');toc
td_model=td(1)+days(tout);
ind1=find(tout>=2,1,'first');
ind2=length(tout);
ind3=find(tout>=7,1,'first');
ind_eval=find(tout>=1,1,'first');

%% Save final states to use as init values for next run
% stateset_init=sewer_out1(end,9:16);
% save stateset_init stateset_init

%% Results
load sewer_out_model

plot_figures(td(ind1:ind2),sewer_out_model(ind1:ind2,7),sewer_out11(ind1:ind2,7),td(ind1),td(ind2),0,7000,'Time (days)','Flow rate (m^3/d)','hydraulic model downstream','conceptual model downstream','conceptualmodel_flowrate_malmo');
plot_figures(td(ind1:ind3),sewer_out_model(ind1:ind3,7),sewer_out11(ind1:ind3,7),td(ind1),td(ind3),0,7000,'Time (days)','Flow rate (m^3/d)','hydraulic model downstream','conceptual model downstream','conceptualmodel_flowrate_snapshot_malmo');

%Temperature
plot_figures3(td(ind1:ind2),T_up_5min.T(ind1:ind2)-273.15,T_down_5min.T(ind1:ind2)-273.15,sewer_out11(ind1:ind2,8)-273.15,td(ind1),td(ind2),10,18,'Time (days)','Temperature ({^\circ}C)','data upstream','data downstream','conceptual model downstream','conceptualmodel_temperature_malmo');
plot_figures3(td(ind1:ind3),T_up_5min.T(ind1:ind3)-273.15,T_down_5min.T(ind1:ind3)-273.15,sewer_out11(ind1:ind3,8)-273.15,td(ind1),td(ind3),10,18,'Time (days)','Temperature ({^\circ}C)','data upstream','data downstream','conceptual model downstream','conceptualmodel_temperature_snapshot_malmo');

%Evaluation metrics
rmse_malmo=rms(T_down_5min.T(ind_eval:ind2)-sewer_out11(ind_eval:ind2,8));
maxerror_malmo= max(abs(T_down_5min.T(ind_eval:ind2)-sewer_out11(ind_eval:ind2,8)));
meanerror_malmo= mean(abs(T_down_5min.T(ind_eval:ind2)-sewer_out11(ind_eval:ind2,8)));
disp(['Root mean squared error conceptual model malmo: ',num2str(rmse_malmo)])
disp(['Max absolute error conceptual model malmo ',num2str(maxerror_malmo)])
disp(['Mean absolute error conceptual model malmo ',num2str(meanerror_malmo)])
