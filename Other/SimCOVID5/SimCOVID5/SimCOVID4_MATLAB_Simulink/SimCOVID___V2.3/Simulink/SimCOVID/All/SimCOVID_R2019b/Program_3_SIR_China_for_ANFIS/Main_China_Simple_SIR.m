% SimCOVID Version 1.0, March 2020
% Author: Ismael Abdulrahman
% This program was made for the paper: "SimCOVID: An Open-Source Simulation Program for the COVID-19 Outbreak"
% The data was taken from: https://www.ecdc.europa.eu/en/geographical-distribution-2019-ncov-cases
% Case: China Outbreak - Simple SIR Model
% Note: Please install the latest version of MATLAB\Simulink (R2019b or newer) 
% as the Simulink programs were built using version R2019b
%==========================================================================
clear all; clc;
N   = 11*10^6; % Population, Susceptible, N, S0
S0  = N;      % Initial value for "Susceptible"
I0  = 1;      % Initial value for "Infectious"
C0  = 1;      % initial value for "Cumulative"
R0  = 0;      % Initial value for "Removed"
E0  = 1;      % Initial value for "Exposed"
D0  = 0;      % Initial value for "Death"
L0  = 2;      % Latent time
st1 = 28;     % starting time 1 (step time)
st2 = 44;     % starting time 2 (step time)
%==========================================================================
%%%%% These are the final values from the optimization obtained by Simulink
gamma = 1.957113604746750; 
beta  = 3.344469065530917;
zeta  = -.034;
I1 = xlsread('COVID-19-geographic-disbtribution-worldwide-2020-04-05','China','E:E');
I  = flip(I1);  
C  =cumsum(I);               % Cumulative sum of I   
t  = 1:length(C);          
xy = [t' C];
xy0 = [t' I];       
%=============================Plotting=====================================
figure (1)
t_length=90;
simOut = sim('Program_3_China_Simple_SIR',t_length);
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
