% SimCOVID Version 1.12, March 2020
% Author: Ismael Abdulrahman
% This program was made for the paper: "SimCOVID: An Open-Source Simulink-Based Program for Simulating the COVID-19 Spread"
% The data was taken from: https://www.ecdc.europa.eu/en/geographical-distribution-2019-ncov-cases
% Case: Italy Outbreak - Variable Gamma, Variable Beta
% The data used for this simulation started from Jan 31 until April 7 (68 days.
%==========================================================================
clc; warning off
P   = 60*10^6; % Population, N
g   = 1;       % Population gain (1 > g > 0)
I0  = 1;       % Initial value for "Infectious"
R0  = 0;       % Initial value for "Removed"
st1 = 40;      % starting time 1 (step time)
st2 = 52;      % starting time 2 (step time)
%====================Loading data and assigning parameters=================
a=0.5;                       % Sigmoid parameter (smooth step function)
gamma1  = 1.80225828092834; 
gamma2  = 0.00100650344260686;
gamma3  = 0.704874058340747;
beta1   = 1.99083773596794;
beta2   = 0.105110391580151;
beta3   = 0.676537806253442;
I1 = xlsread('COVID-19-geographic-disbtribution-worldwide-2020-04-07','Italy','E:E');
Im  = flip(I1);  Im=Im(32:end); %  The first 31 elements are zeros
Cm  = cumsum(Im);               % Cumulative sum of I   
t   = 1:length(Cm);          
xy  = [t' Cm];
xy0 = [t' Im];
% st3 = length(Im);
st3 = 0;
%%%%===========================Prepare Parameter Inputs=======================
t_length=180;                       % number of days
isModelOpen = bdIsLoaded('SIR_Italy_VariableBetaGamma100'); % We need this later on to close the parallel loop
open_system('SIR_Italy_VariableBetaGamma100');
beta_new=1;
beta_step = beta_new*[0.995 1 1.005];   % These are normalized beta changes (10%)
numSims = length(beta_step);
for i = numSims:-1:1
    in(i) = Simulink.SimulationInput('SIR_Italy_VariableBetaGamma100');
    in(i) = in(i).setBlockParameter('SIR_Italy_VariableBetaGamma100/Value2', 'Value', num2str(beta_step(i)));
end
in = in.setBlockParameter('SIR_Italy_VariableBetaGamma100/Sigmoidal MF1','a',num2str(-a));
in = in.setBlockParameter('SIR_Italy_VariableBetaGamma100/Sigmoidal MF1','c',num2str(st1));
in = in.setBlockParameter('SIR_Italy_VariableBetaGamma100/Sigmoidal MF2','a',num2str(a));
in = in.setBlockParameter('SIR_Italy_VariableBetaGamma100/Sigmoidal MF2','c',num2str(st1));
in = in.setBlockParameter('SIR_Italy_VariableBetaGamma100/Sigmoidal MF3','a',num2str(a));
in = in.setBlockParameter('SIR_Italy_VariableBetaGamma100/Sigmoidal MF3','c',num2str(st2));
in = in.setBlockParameter('SIR_Italy_VariableBetaGamma100/gamma1','Gain',num2str(gamma1));
in = in.setBlockParameter('SIR_Italy_VariableBetaGamma100/gamma2','Gain',num2str(gamma2));
in = in.setBlockParameter('SIR_Italy_VariableBetaGamma100/gamma3','Gain',num2str(gamma3));
in = in.setBlockParameter('SIR_Italy_VariableBetaGamma100/beta1','Gain',num2str(beta1));
in = in.setBlockParameter('SIR_Italy_VariableBetaGamma100/beta2','Gain',num2str(beta2));
in = in.setBlockParameter('SIR_Italy_VariableBetaGamma100/beta3','Gain',num2str(beta3));
in = in.setBlockParameter('SIR_Italy_VariableBetaGamma100/g','Gain',num2str(g));
in = in.setBlockParameter('SIR_Italy_VariableBetaGamma100/Value1','Value',num2str(1));
in = in.setBlockParameter('SIR_Italy_VariableBetaGamma100/Switch','Threshold',num2str(st3));
in = in.setBlockParameter('SIR_Italy_VariableBetaGamma100/Population','Value',num2str(P));
in = in.setBlockParameter('SIR_Italy_VariableBetaGamma100/I0','Value',num2str(I0));
in = in.setBlockParameter('SIR_Italy_VariableBetaGamma100/R0','Value',num2str(R0));
% Use out(1).ErrorMessage for any parallel computing problem
%=========================Running Parallel Loop============================
out = parsim(in,'ShowProgress','on'); % 'ShowSimulationManager','on',
t=out(1).tout;
I1=out(1).I; I2=out(2).I; I3=out(3).I;
C1=out(1).C; C2=out(2).C; C3=out(3).C;
tt = vertcat(t,flipud(t));
II = vertcat(I1,flipud(I3));
CC = vertcat(C1,flipud(C3));
%==============================Plotting C==================================
figure('Color',[1 1 1],'units','normalized','outerposition',[0 0 1 1])
plot(t,C1, 'Color','[1 1 1]'); hold on
plot(t,C3, 'Color','[1 1 1]')
patch(tt,CC,[0.75 0.75 0.75],'FaceAlpha',0.45, 'EdgeColor' , [1 1 1]) % trancperency (0 is for complete transp and 1 for fully opaque
hold on
plot(t,C2,'b','LineWidth',3,'Color', '[0, 0.4470, 0.7410]')
hold on
plot(Cm,'*','LineWidth',1,'LineStyle',':','MarkerSize',10, 'Color','[0.8500, 0.3250, 0.0980]')
% Use the color codes from this site: http://math.loyola.edu/~loberbro/matlab/html/colorsInMatlab.html
grid on;set(gca,'fontsize',16);
xlabel('Time (day)','FontSize',16,'FontWeight','bold');
xlim([0 t_length])
ylabel('Cumulative Infectious (person)','FontSize',16,'FontWeight','bold');
dateaxis('x', 6, '31-Jan-2020')
%==============================Plotting I==================================
figure('Color',[1 1 1],'units','normalized','outerposition',[0 0 1 1])
plot(t,I1, 'Color','[1 1 1]'); hold on
plot(t,I3, 'Color','[1 1 1]')
patch(tt,II,[0.75 0.75 0.75],'FaceAlpha',0.45, 'EdgeColor' , [1 1 1]) % trancperency (0 is for complete transp and 1 for fully opaque
hold on
plot(t,I2,'b','LineWidth',3,'Color', '[0, 0.4470, 0.7410]')
hold on
plot(Im,'*','LineWidth',1,'LineStyle',':','MarkerSize',10, 'Color','[0.8500, 0.3250, 0.0980]')
% Use the color codes from this site: http://math.loyola.edu/~loberbro/matlab/html/colorsInMatlab.html
grid on;set(gca,'fontsize',16);
xlabel('Time (day)','FontSize',16,'FontWeight','bold');
xlim([0 t_length])
ylabel('Daily Infectious (person)','FontSize',16,'FontWeight','bold');
dateaxis('x', 6, '31-Jan-2020')
%===========================Close the parallel loop========================
if(~isModelOpen)
    close_system('SIR_Italy_VariableBetaGamma100', 0);
end
poolobj = gcp('nocreate');
% % ==========================================================================
