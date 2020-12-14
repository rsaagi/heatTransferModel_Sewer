
function [hAxes] = plot_figures(x,y1,y2,xlim0,xlim1,ylim0,ylim1,label1,label2,legend1,legend2,fig_name)
%plot function for two y variables and one x variable. The plot is set to a
%specific size and set of colors for easy use in publications
hFig = figure;

hAxes = axes('Parent',hFig,...
    'FontName','Arial','FontSize',9,'XGrid','off', 'YGrid','off', 'Box', 'on');
hold('on');

plot(hAxes, x,y1,'LineStyle', '-', 'Color', [0/256 130/256 200/256],... %blueish
    'LineWidth', 1);
plot(hAxes, x,y2,'LineStyle', '-', 'Color', [230/256 25/256 75/256],... %reddish
    'LineWidth', 1);



xlabel(label1,'FontName','Arial');
ylabel(label2,'FontName','Arial');
yticks([linspace(ylim0,ylim1,5)]);
datetick('x','dd-mmm','Keepticks');
hLegend=legend(hAxes, legend1, legend2);

set(hLegend, 'FontSize',9,'Box','off','FontName','Arial');
set(hAxes,'XLim',[xlim0 xlim1],'YLim',[ylim0,ylim1]);

fig = gcf;
fig.PaperUnits = 'centimeters';
fig.PaperPosition = [0 0 9 0.7*9];

print(fig_name,'-dmeta')