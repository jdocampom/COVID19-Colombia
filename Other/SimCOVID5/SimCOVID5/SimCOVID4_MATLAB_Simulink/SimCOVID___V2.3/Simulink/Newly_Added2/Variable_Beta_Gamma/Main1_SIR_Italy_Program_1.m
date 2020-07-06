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
% %========================Plotting basic I and C==========================
t_length = 180;
simOut   = sim('SIR_Italy_VariableBetaGamma100',t_length);
Cs       = simOut.C;
Is       = simOut.I;
ts       = simOut.tout;

figure('Color',[1 1 1],'units','normalized','outerposition',[0 0 1 1])
plot(ts,Cs, 'linewidth',3)
grid on;set(gca,'fontsize',16);
xlabel('Time (day)','FontSize',16,'FontWeight','bold');
xlim([0 t_length])
ylabel('Cumulative Infectious (person)','FontSize',16,'FontWeight','bold');
hold on
plot(Cm,'*','LineWidth',1,'LineStyle',':','MarkerSize',10, 'Color','[0.8500, 0.3250, 0.0980]')
legend({'simulated infectious','reported cases'}, 'FontSize',12);
dateaxis('x', 6, '31-Jan-2020')

figure('Color',[1 1 1],'units','normalized','outerposition',[0 0 1 1])
plot(ts,Is, 'linewidth',3)
grid on;set(gca,'fontsize',16);
xlabel('Time (day)','FontSize',16,'FontWeight','bold');
xlim([0 t_length])
ylabel('Cumulative Infectious (person)','FontSize',16,'FontWeight','bold');
hold on
plot(Im,'*','LineWidth',1,'LineStyle',':','MarkerSize',10, 'Color','[0.8500, 0.3250, 0.0980]')
legend({'simulated infectious','reported cases'}, 'FontSize',12);
dateaxis('x', 6, '31-Jan-2020')
%%%==========================================================================
