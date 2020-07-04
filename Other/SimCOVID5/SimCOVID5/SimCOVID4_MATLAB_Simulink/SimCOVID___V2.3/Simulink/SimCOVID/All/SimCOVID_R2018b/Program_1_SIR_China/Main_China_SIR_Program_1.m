% SimCOVID Version 1.0, March 2020
% Author: Ismael Abdulrahman
% This program was made for the paper: "SimCOVID: An Open-Source Simulation Program for the COVID-19 Outbreak"
% The data was taken from: https://www.ecdc.europa.eu/en/geographical-distribution-2019-ncov-cases
% Case: China Outbreak
% Note: Please install the latest version of MATLAB\Simulink (R2019b or newer) 
% as the Simulink programs were built using version R2019b
%==========================================================================
clear all; clc;
P   = 11*10^6; % Population, Susceptible, N, S0
I0  = 1;       % Initial value for "Infectious"
R0  = 0;       % Initial value for "Removed"
st1 = 28;      % starting time 1 (step time)
st2 = 44;      % starting time 2 (step time)
a   = 0.5;     % parameter of the sigmoid function
%==========================================================================
%%%%% These are the final values from the optimization obtained by Simulink
beta    = 1.17386114378623;
gamma1  = 0.65905784176993; 
gamma2  = 0.922873642259164;
gamma3  = 1.06040119904744;
zeta1   = -0.251586911587193;
zeta2   = -0.124211246181375;
zeta3   = -0.311839590070327;
I1  = xlsread('COVID-19-geographic-disbtribution-worldwide-2020-04-05','China','E:E');
I   = flip(I1);  
C   = cumsum(I);               % Cumulative sum of I   
t   = 1:length(C);          
xy  = [t' C];
xy0 = [t' I];
%=============================Plotting=====================================
figure (1)
t_length=90;
simOut = sim('Program_1_SIR_China_R2018b',t_length);
plot(simOut.C, 'linewidth',2)
grid on;set(gca,'fontsize',16);
xlabel('Time (day)','FontSize',16,'FontWeight','bold');
xlim([0 t_length])
ylabel('Cumulative Infectious (person)','FontSize',16,'FontWeight','bold');
hold on
plot(C, 'linewidth',2)
legend({'simulated infectious','reported cases'}, 'FontSize',12);
dateaxis('x', 6, '31-Dec-2020')

figure (2)
plot(simOut.I, 'linewidth',2)
grid on;set(gca,'fontsize',16);
xlabel('Time (day)','FontSize',16,'FontWeight','bold');
xlim([0 t_length])
ylabel('Cumulative Infectious (person)','FontSize',16,'FontWeight','bold');
hold on
plot(I, 'linewidth',2)
legend({'simulated infectious','reported cases'}, 'FontSize',12);
dateaxis('x', 6, '31-Dec-2020')
%==========================================================================
