% SimCOVID Version 1.0, March 2020
% Author: Ismael Abdulrahman
% This program was made for the paper: "SimCOVID: An Open-Source Simulation Program for the COVID-19 Outbreak"
% The data was taken from: https://www.ecdc.europa.eu/en/geographical-distribution-2019-ncov-cases
% Case: China Outbreak - ANFIS
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
%==========================================================================