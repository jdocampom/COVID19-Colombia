% SimCOVID Version 1.0, March 2020
% Author: Ismael Abdulrahman
% This program was made for the paper: "SimCOVID: An Open-Source Simulation Program for the COVID-19 Outbreak"
% The data was taken from: https://www.ecdc.europa.eu/en/geographical-distribution-2019-ncov-cases
% Case: Italy Outbreak - SEIRD
% Note: Please install the latest version of MATLAB\Simulink (R2019b or newer) 
% as the Simulink programs were built using version R2019b
%==========================================================================
clear all; clc;
N   = 60*10^6; % Population, Susceptible, N, S0
S0  = N;      % Initial value for "Susceptible"
I0  = 1;      % Initial value for "Infectious"
R0  = 0;      % Initial value for "Removed"
E0  = 1;      % Initial value for "Exposed"
D0  = 1;      % Initial value for "Death"
L0  = 2;      % Latent time
st1 = 40;     % starting time 1 (step time)
st2 = 52;     % starting time 2 (step time)
%==========================================================================
%%%%% These are the final values from the optimization obtained by Simulink
beta   = 1.971574662033975;
gamma1 = 0.019992964106553; 
gamma2 = 0.015981365215522;
gamma3 = 0.564739764641216;
zeta1  = -2.179313856128234;
zeta2  = -2.080060703203781;
zeta3  = -1.046051249951776;
I = xlsread('Italy','Italy2','B:B');
C = xlsread('Italy','Italy2','D:D');
% Note: to calculate the D/I ratio in Simulink, I-variable must be a nonzero variable. Therefore, some
% initial days with zero infectious are changed wto one
t  = 1:length(C);          
xy = [t' C];
xy0 = [t' I];
%=============================Plotting=====================================
figure (1)
t_length = 70;
simOut = sim('Program_2_Italy_SEIRD_R2019a',t_length);
plot(simOut.C, 'linewidth',2)
grid on;set(gca,'fontsize',16);
xlabel('Time (day)','FontSize',16,'FontWeight','bold');
xlim([0 t_length])
ylabel('Cumulative Infectious (person)','FontSize',16,'FontWeight','bold');
hold on
plot(C, 'linewidth',2)
legend({'simulated infectious','reported cases'}, 'FontSize',12);
dateaxis('x', 6, '31-Jan-2020')

figure (2)
plot(simOut.I, 'linewidth',2)
grid on;set(gca,'fontsize',16);
xlabel('Time (day)','FontSize',16,'FontWeight','bold');
xlim([0 t_length])
ylabel('Infectious (person)','FontSize',16,'FontWeight','bold');
hold on
plot(I, 'linewidth',2)
legend({'simulated infectious','reported cases'}, 'FontSize',12);
dateaxis('x', 6, '31-Jan-2020')
%==========================================================================
