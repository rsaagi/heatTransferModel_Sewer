%% Mechanistic model
run_sewermodel_malmo_mechanistic
sewer_out_model=sewer_out11;
save('sewer_out_model.mat','sewer_out_model')
clearvars

%% Conceptual model
run_sewermodel_malmo_conceptual

%% Additional figures comparing mechanistic and conceptual model

%Temp mechanistic vs conceptual
plot_figures3(td(ind1:ind2),T_down_5min.T(ind1:ind2)-273.15,sewer_out_model(ind1:ind2,8)-273.15,sewer_out11(ind1:ind2,8)-273.15,td(ind1),td(ind2),10,18,'Time (days)','Temperature ({^\circ}C)','data downstream','mechanistic model downstream','conceptual model downstream','modelcomparison_temperature_malmo');
plot_figures3(td(ind1:ind3),T_down_5min.T(ind1:ind3)-273.15,sewer_out_model(ind1:ind3,8)-273.15,sewer_out11(ind1:ind3,8)-273.15,td(ind1),td(ind3),10,18,'Time (days)','Temperature ({^\circ}C)','data downstream','mechanistic model downstream','conceptual model downstream','modelcomparison_temperature_snapshot_malmo');