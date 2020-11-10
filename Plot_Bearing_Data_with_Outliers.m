%% 
clc
clear all
close all
addpath(genpath(fileparts(mfilename('fullpath'))));
rng('default')
rng(1)

%% generate the simulation signal for RobustTrend
load RMS1.mat
x = Full_RMS(1:2:end);
x = x(:);
y = x;
y(500) = y(499) + 0.2;
y(600) = y(600) + 0.18;
n = length(x);
points = (0:n-1) * 20 / 60;

%%%%%%%%% robust
k = 2;
lambda = 990;
[x_robust] = Robust_ETF(y, lambda, k);


%% Figure initialization
global FontSize FontName;
Tstring =  'Times (h)'; 
Fstring = 'Frequency (Hz)';
Astring =  'RMS';
FontSize = 9;   FontName = 'Times New Roman';
MarkerSize = 4;  LineWidth = 1;
%%
FlagFigureAutoSave = 1;
currentFolder = pwd;



%% Print the time domain
[WindowPosition,h1] = Subfigure11_cm(8.6, 2, 1.3, 0.3, 0.2, 0.8);
figure(4);clf;
set(gcf, 'NumberTitle','off','Name','7d_Robust_Experiment');
set(gcf, 'Units', 'centimeters');
set(gcf,'position',WindowPosition);
set(gcf, 'PaperPositionMode', 'auto');   
h1_ap=axes('position',h1); 
hold on
ph(1) = plot(points, x, 'b-', 'LineWidth', LineWidth);
ph(2) = plot(points, x_robust, 'r-', 'LineWidth', LineWidth);
hold off
legend1 = legend(ph, 'Original', 'RobustETF');
set(legend1,'location','southeast','Orientation','horizontal', 'FontSize',FontSize,'FontName',FontName, 'Interpreter', 'latex')
legend boxoff

xlim_min = 0; xlim_max = max(points);
ylim_min = min(x);            ylim_max = max(abs(x));  
xylim = [xlim_min,xlim_max,ylim_min,ylim_max]; axis(xylim);

LabelX_Linchao(h1_ap,Tstring,xylim,0.22);
LabelY_Linchao(h1_ap,Astring,xylim,0.1);
set(h1_ap,'FontSize',FontSize,'FontName',FontName);

annotation('textbox',[0 1 0.03 0.03],'String',{'(d)'},'FontSize',FontSize+1,'FontName',FontName,'FontWeight','bold','FitBoxToText','off','LineStyle','none');
% save figure
SaveFigureLinchao('Bearing_Data_with_Outliers',FlagFigureAutoSave,currentFolder)



