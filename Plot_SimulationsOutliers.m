%% 
clc
clear all
close all
addpath(genpath(fileparts(mfilename('fullpath'))));
rng('default')
rng(1)

%% Figure initialization
global FontSize FontName;
Tstring =  'Points'; 
Fstring = 'Frequency (Hz)';
Astring =  'Amp';
FontSize = 9;   FontName = 'Times New Roman';
MarkerSize = 4;  LineWidth = 1;
%%
FlagFigureAutoSave = 1;
currentFolder = pwd;
%% generate the simulation signal for RobustTrend
x_first = zeros(200, 1);
x_first(21:40) = 1.5;
x_first(41:70) = -1.5;
x_first(71:100) = 1;
x_first(101:160) = 0;
x_first(161:170) = 1;
x_first(171:180) = -1;
x_first(181:200) = 0;
x_first = x_first - mean(x_first);
%sin_wave = sin(2*pi*(1:200)/200);
x_second = interp1([1 50 75 130 150 200], 2 * [0.3 1 -0.2 -0.5 0.9 0.6],1:200)';
x_second = x_second - mean(x_second);

points = 1 : 200;
n = length(points);
sigma = 0.06;
rand_index = randperm(n);
rand_outlier = rand_index(1 : round(n*sigma));
y_second = x_second + 0.2 * randn(n, 1);
y_second(rand_outlier) = y_second(rand_outlier) + 2 * randsrc(round(n*sigma), 1);

y_first = x_first + 0.2 * randn(n, 1);
y_first(rand_outlier) = y_first(rand_outlier) + 2 * randsrc(round(n*sigma), 1);

SNR_first = C_SNR( x_first , y_first )

SNR_first = C_SNR( x_second , y_second )
%% Print the time domain
[WindowPosition,h1] = Subfigure11_cm(8.6, 2, 1, 0.3, 0.2, 0.8);
figure(1);clf;
set(gcf, 'NumberTitle','off','Name','5aSimulation_First');
set(gcf, 'Units', 'centimeters');
set(gcf,'position',WindowPosition);
set(gcf, 'PaperPositionMode', 'auto');   
h1_ap=axes('position',h1); 
hold on
ph(1) = plot(points, y_first, 'b-', 'LineWidth', LineWidth);
ph(2) = plot(points, x_first, 'r--', 'LineWidth', LineWidth);
hold off
% box on
legend1 = legend(ph, 'Original', 'Trend');
set(legend1,'location','best','Orientation','horizontal', 'FontSize',FontSize,'FontName',FontName, 'Interpreter', 'latex')
legend boxoff

xlim_min = 0; xlim_max = max(points);
ylim_min = -max(abs(y_first))*1.1;            ylim_max = max(abs(y_first))*1.3;  
xylim = [xlim_min,xlim_max,ylim_min,ylim_max]; axis(xylim);

LabelX_Linchao(h1_ap,Tstring,xylim,0.22);
LabelY_Linchao(h1_ap,Astring,xylim,0.05);
set(h1_ap,'FontSize',FontSize,'FontName',FontName);

annotation('textbox',[0 1 0.03 0.03],'String',{'(a)'},'FontSize',FontSize+1,'FontName',FontName,'FontWeight','bold','FitBoxToText','off','LineStyle','none');
% save figure
SaveFigureLinchao('Simulation_First',FlagFigureAutoSave,currentFolder)
%% Print the time domain

figure(2);clf;
set(gcf, 'NumberTitle','off','Name','5f_Simulation_Second');
set(gcf, 'Units', 'centimeters');
set(gcf,'position',WindowPosition);
set(gcf, 'PaperPositionMode', 'auto');   
h1_ap=axes('position',h1); 
hold on
ph(1) = plot(points, y_second, 'b-', 'LineWidth', LineWidth);
ph(2) = plot(points, x_second, 'r--', 'LineWidth', LineWidth);
hold off
legend1 = legend(ph, 'Original', 'Trend');
set(legend1,'location','best','Orientation','horizontal', 'FontSize',FontSize,'FontName',FontName, 'Interpreter', 'latex')
legend boxoff

xlim_min = 0; xlim_max = max(points);
ylim_min = min(y_second)*1.1;            ylim_max = max(abs(y_second))*1.2;  
xylim = [xlim_min,xlim_max,ylim_min,ylim_max]; axis(xylim);

LabelX_Linchao(h1_ap,Tstring,xylim,0.22);
LabelY_Linchao(h1_ap,Astring,xylim,0.05);
set(h1_ap,'FontSize',FontSize,'FontName',FontName);

annotation('textbox',[0 1 0.03 0.03],'String',{'(a)'},'FontSize',FontSize+1,'FontName',FontName,'FontWeight','bold','FitBoxToText','off','LineStyle','none');
% save figure
SaveFigureLinchao('Simulation_Second',FlagFigureAutoSave,currentFolder)



%% First
k = 1;
load First_Outliers.mat
%%%%%%%%% robust
[~, index] = min(MAE_robust(4, :, 1));
lambda = lam(index);
[x_robust] = Robust_ETF(y_first, lambda, k);
[MSE_robust, MAE_robust] = Indicator(x_first, x_robust);

%% Print the time domain
[WindowPosition,h1] = Subfigure11_cm(8.6, 2, 1, 0.3, 0.4, 0.8);
figure(6);clf;
set(gcf, 'NumberTitle','off','Name','5e_Robust_First');
set(gcf, 'Units', 'centimeters');
set(gcf,'position',WindowPosition);
set(gcf, 'PaperPositionMode', 'auto');   
h1_ap=axes('position',h1); 
hold on
ph(1) = plot(points, x_robust, 'b-', 'LineWidth', LineWidth);
ph(2) = plot(points, x_first, 'r--', 'LineWidth', LineWidth);
hold off
legend1 = legend(ph, 'RobustETF', 'Trend');
set(legend1,'location','best','Orientation','horizontal', 'FontSize',FontSize,'FontName',FontName, 'Interpreter', 'latex')
legend boxoff
title(['MSE=', num2str(round(MSE_robust*10000)/10000), ', MAE=', num2str(round(MAE_robust*10000)/10000)],'FontSize',FontSize,'FontName',FontName);

xlim_min = 0; xlim_max = max(points);
ylim_min = min(x_robust)*1.1;            ylim_max = max(abs(x_robust))*1.5;  
xylim = [xlim_min,xlim_max,ylim_min,ylim_max]; axis(xylim);

LabelX_Linchao(h1_ap,Tstring,xylim,0.22);
LabelY_Linchao(h1_ap,Astring,xylim,0.05);
set(h1_ap,'FontSize',FontSize,'FontName',FontName);

annotation('textbox',[0 1 0.03 0.03],'String',{'(e)'},'FontSize',FontSize+1,'FontName',FontName,'FontWeight','bold','FitBoxToText','off','LineStyle','none');
% save figure
SaveFigureLinchao('Simulation_Robust_First',FlagFigureAutoSave,currentFolder)



%% Second
k = 2;
load Second_Outliers.mat
%%%%%%%%% robust
[~, index] = min(MAE_robust(4, :, 1));
lambda = lam(index);
[x_robust] = Robust_ETF(y_second, lambda, k);
[MSE_robust, MAE_robust] = Indicator(x_second, x_robust);

%% Print the time domain
[WindowPosition,h1] = Subfigure11_cm(8.6, 2, 1, 0.3, 0.4, 0.8);
figure(10);clf;
set(gcf, 'NumberTitle','off','Name','5j_Robust_Second');
set(gcf, 'Units', 'centimeters');
set(gcf,'position',WindowPosition);
set(gcf, 'PaperPositionMode', 'auto');   
h1_ap=axes('position',h1); 
hold on
ph(1) = plot(points, x_robust, 'b-', 'LineWidth', LineWidth);
ph(2) = plot(points, x_second, 'r--', 'LineWidth', LineWidth);
hold off
legend1 = legend(ph, 'RobustETF', 'Trend');
set(legend1,'location','best','Orientation','horizontal', 'FontSize',FontSize,'FontName',FontName, 'Interpreter', 'latex')
legend boxoff
title(['MSE=', num2str(round(MSE_robust*10000)/10000), ', MAE=', num2str(round(MAE_robust*10000)/10000)],'FontSize',FontSize,'FontName',FontName);

xlim_min = 0; xlim_max = max(points);
ylim_min = min(x_robust)*1.1;            ylim_max = max(abs(x_robust))*1.5;  
xylim = [xlim_min,xlim_max,ylim_min,ylim_max]; axis(xylim);

LabelX_Linchao(h1_ap,Tstring,xylim,0.22);
LabelY_Linchao(h1_ap,Astring,xylim,0.05);
set(h1_ap,'FontSize',FontSize,'FontName',FontName);

annotation('textbox',[0 1 0.03 0.03],'String',{'(e)'},'FontSize',FontSize+1,'FontName',FontName,'FontWeight','bold','FitBoxToText','off','LineStyle','none');
% save figure
SaveFigureLinchao('Simulation_Robust_Second',FlagFigureAutoSave,currentFolder)


