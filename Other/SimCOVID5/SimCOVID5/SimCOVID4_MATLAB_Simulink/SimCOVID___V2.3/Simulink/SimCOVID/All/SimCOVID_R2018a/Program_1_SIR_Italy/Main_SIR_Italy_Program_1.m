% SimCOVID Version 1.0, March 2020
% Author: Ismael Abdulrahman
% This program was made for the paper: "SimCOVID: An Open-Source Simulation Program for the COVID-19 Outbreak"
% The data was taken from: https://www.ecdc.europa.eu/en/geographical-distribution-2019-ncov-cases
% Case: Italy Outbreak
% Note: Please install the latest version of MATLAB\Simulink (R2019b or newer) 
% as the Simulink programs were built using version R2019break
%==========================================================================
clear all; clc;
P   = 60*10^6; % Population, Susceptible, N, S0
I0  = 1;       % Initial value for "Infectious"
R0  = 0;       % Initial value for "Removed"
st1 = 40;      % starting time 1 (step time)
st2 = 52;      % starting time 2 (step time)
%==========================================================================
a=0.5;
gamma=0.05;
beta    = 1.53602767817435;
gamma1  =  0.0128250346772537; 
gamma2  = 0.00483679113978549;
gamma3  = 0.149311333504727;
zeta1   = -2.03011707682403;
zeta2   = -2.60970797840359;
zeta3   = -2.55914577598691;
I1 = xlsread('COVID-19-geographic-disbtribution-worldwide-2020-04-07','Italy','E:E');
Im  = flip(I1);  Im=Im(32:end); %  The first 31 elements are zeros
Cm  =cumsum(Im);               % Cumulative sum of I   
t  = 1:length(Cm);          
xy = [t' Cm];
xy0 = [t' Im];
%=============================Plotting=====================================
t_length=180;
simOut = sim('Program_1_SIR_Italy_R2018a',t_length);
Cs=simOut.C;
Is=simOut.I;
ts=simOut.tout;

figure (1)
plot(ts,Cs, 'linewidth',2)
grid on;set(gca,'fontsize',16);
xlabel('Time (day)','FontSize',16,'FontWeight','bold');
xlim([0 t_length])
ylabel('Cumulative Infectious (person)','FontSize',16,'FontWeight','bold');
hold on
plot(Cm, 'linewidth',2)
legend({'simulated infectious','reported cases'}, 'FontSize',12);
dateaxis('x', 6, '31-Jan-2020')

figure (2)
plot(ts,Is, 'linewidth',2)
grid on;set(gca,'fontsize',16);
xlabel('Time (day)','FontSize',16,'FontWeight','bold');
xlim([0 t_length])
ylabel('Cumulative Infectious (person)','FontSize',16,'FontWeight','bold');
hold on
plot(Im, 'linewidth',2)
legend({'simulated infectious','reported cases'}, 'FontSize',12);
dateaxis('x', 6, '31-Jan-2020')
%%%%==========================================================================

