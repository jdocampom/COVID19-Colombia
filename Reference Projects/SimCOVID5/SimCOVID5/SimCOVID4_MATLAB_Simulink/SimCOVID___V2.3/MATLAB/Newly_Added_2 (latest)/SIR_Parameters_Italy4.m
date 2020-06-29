% SimCOVID Version 2.0, MATLAB Part, April 25, 2020
% Author: Ismael Abdulrahman
% The data was taken from (updated daily): https://www.ecdc.europa.eu/en/geographical-distribution-2019-ncov-cases
% Case: Italy Outbreak
% Note 1: You can control your curve-fitting by changing the values of st1, st2, and a. Choose critical points for st1 and st2 such as inflection points of the curve (Cm or Im). If not, choose st1 and st2 at one-third and two-third time of simulation (try it lastly). You can also add or remove some of the beta function branches. Choose "a" between 0 and 1 as initial values, if possible. Low values of "a" give a smoother curve.
% Note 2: The selected initial values are so so important. You can help the optimizer tool to find the solution much much faster if you choose these values wisely.
% Note 3: The lower bound for the parameters here must give positive beta and gamma. Changing the upper bound limits can affect the solution and the time to solve the problem.
% Note 4: This program uses SIR model. However, for the cumulative infectious, an integrator is added to cumulatively sum I, and the total number of DEs becomes 4. If you want to use "cumsum" to integrate I with 2 DEs, use it without customized step size (see the other files).
% Note 5: This case is for a variable beta, variable gamma with 9 parameters to be estimated
function SIR_Parameters_Italy4
% function [Parameters] = SIR_Parameters_Italy4
I0  = 3;                                    % Inital value of I (number of infectious befor simulation start time
C0  = 3;                                    % Initial value of C (cumulative infectious)
R0  = 0;                                    % Initial value of R "Recovered"
N   = 60*10^6;                              % Population
S0  = N-I0;                                 % Initial value of "Susceptible"
% st1 = 50; st2 = 55;                       % Step times for the sigmoid function (chosen empirically by observing Cm plot)
x0  = [S0; I0; R0; C0 ];                    % Initial values for the states (S,I,R) used in 'ode45'                  % To give a smooth curve, choose a step size less than 1 (0.1 is preferred)
Simulation_length = 180;                    % Simulation time (day)
step_size = 0.1;                            % To give a smooth curve, choose a step size less than 1 (0.1 is preferred)
Simulated_t = 1:step_size:Simulation_length;
betag = [0.98 1.0 1.02];                  % Beta change at time zero.
Im0   = xlsread('COVID-19-geographic-disbtribution-worldwide-2020-04-25','Italy','E:E');
Im    = flip(Im0);  Im=Im(32:end);          %  The first 31 elements are zeros
Cm    = cumsum(Im);                         % Cumulative sum of I
tm    = 1:length(Cm);                       % Time (day)
[Parameters] = SIR_Model;                   % Run the function "SIR_Model" to get the basic estimated parameters without change in the beta function
%==========================Solving SIR Three Times=========================
Simulated_SIR=zeros(length(Simulated_t),3);
for i = 1:3
    [tout,Simulated_SIR(:,i)] = SolveSIR2(Parameters, Simulated_t); % Solve SIR with three values of beta (standard, plus minus a percentage)
end
%==============================Plotting I's================================
t=tout;
I1 = Simulated_SIR(:,1); I2 = Simulated_SIR(:,2); I3 = Simulated_SIR(:,3);
tt = vertcat(t,flipud(t));
II = vertcat(I1,flipud(I3));
figure('Color',[1 1 1],'units','normalized','outerposition',[0 0 1 1])
plot(t,I1, 'Color','[1 1 1]'); hold on
plot(t,I3, 'Color','[1 1 1]')
patch(tt,II,[0.75 0.75 0.75],'FaceAlpha',0.45, 'EdgeColor' , [1 1 1]) % trancperency (0 is for complete transp and 1 for fully opaque
hold on
plot(t,I2,'b','LineWidth',3,'Color', '[0, 0.5, 0]')
hold on
plot(Im,'p','LineWidth',1,'LineStyle',':','MarkerSize',10, 'Color','[0, 0.4470, 0.7410]')
% Use the color codes from this site: http://math.loyola.edu/~loberbro/matlab/html/colorsInMatlab.html
grid on;set(gca,'fontsize',16); grid minor;
xticks(0:Simulation_length/10:Simulation_length); 
xlabel('Time (day)','FontSize',16,'FontWeight','bold');
xlim([0 Simulation_length])
ylabel('Daily Infectious (person)','FontSize',16,'FontWeight','bold');
dateaxis('x', 6, '31-Jan-2020')
%===========================Printing SIR Parameters========================
fprintf(1,'\tRate Constants:\n')
for i = 1:length(Parameters)
    fprintf(1, '\t\tParameters(%d) = %8.5f\n', i, Parameters(i))
end
%======================Curve Fitting Function (Called)=====================
    function [Parameters] = SIR_Model
        Parameters0 = [0.5; 0.4; 0.2; 0; 50; 0.5; 0.2; 0.2; 0.2];      % Initial values for 'lsqcurvefit'
        options = optimoptions('lsqcurvefit','Algorithm','trust-region-reflective','Display', 'iter','FunctionTolerance',1e-15,'MaxFunctionEvaluations',100000*2,'MaxIterations',10000,'OptimalityTolerance',1e-15); % levenberg-marquardt
        lb = [ 0; 0; 0; 0; 0; 0; 0; 0; 0];                             % minimum values for the parameters (a must)
        ub = 100*[ 1; 1; 1; 1; 1; 0.01; 1; 1; 1];                      % maximum values for the parameters (optional)
        [Parameters] = lsqcurvefit(@SolveSIR1,Parameters0,tm,Cm, lb, ub, options);    % Least-Square curve fitting
    end
%============================= Sigmoid Function============================
    function s = sigmoid(t,c,a)
        s = 1./(1 + exp(-a.*(t-c)));
    end
%============Solving SIR with One Output (Needed for lsqcurvefit)==========
    function SIRout1 = SolveSIR1(Parameters,t)
        [~,SIRout0] = ode45(@DEs,t,x0);
        function dSIR=DEs(t,x)
            betaf1  = Parameters(1)*sigmoid(t,Parameters(4),-Parameters(6));   % See the Simulink file for this function
            betaf2  = Parameters(2)*abs(sigmoid(t,Parameters(4),Parameters(6)) - sigmoid(t,Parameters(5),Parameters(6)));
            betaf3  = Parameters(3)*sigmoid(t,Parameters(5),Parameters(6));
            betaf   = betaf1+betaf2+betaf3;
            gammaf1  = Parameters(7)*sigmoid(t,Parameters(4),-Parameters(6));   % See the Simulink file for this function
            gammaf2  = Parameters(8)*abs(sigmoid(t,Parameters(4),Parameters(6)) - sigmoid(t,Parameters(5),Parameters(6)));
            gammaf3  = Parameters(9)*sigmoid(t,Parameters(5),Parameters(6));
            gammaf   = gammaf1+gammaf2+gammaf3;
            dxdt    = zeros(3,1);                                    % x stands for "state" not x-axis, t is on x-axis
            dxdt(1) = -betaf/N*x(1)*x(2);                            % DE-1: S  (beta has now a new gain "betag")
            dxdt(2) = betaf/N*x(1)*x(2)-gammaf*x(2);                 % DE-2: I  (beta has now a new gain "betag")
            dxdt(3) = gammaf*x(2);                                   % DE-3: R
            dxdt(4) = x(2);                                          % DE-4: C (you can eliminate this state by using "cumsum(I)' but for different step sizes, it might not be the best option)
            dSIR    = dxdt;
        end
        SIRout1 = SIRout0(:,4) ;
    end
%=================Solving SIR with All Outputs (for Plotting)==============
    function [tout, SIRout1] = SolveSIR2(Parameters,t)
        [tout,SIRout0] = ode45(@DEs,t,x0);
        function dSIR=DEs(t,x)
            betaf1  = Parameters(1)*sigmoid(t,Parameters(4),-Parameters(6));   % See the Simulink file for this function
            betaf2  = Parameters(2)*abs(sigmoid(t,Parameters(4),Parameters(6)) - sigmoid(t,Parameters(5),Parameters(6)));
            betaf3  = Parameters(3)*sigmoid(t,Parameters(5),Parameters(6));
            betaf   = betaf1+betaf2+betaf3;
            gammaf1  = Parameters(7)*sigmoid(t,Parameters(4),-Parameters(6));   % See the Simulink file for this function
            gammaf2  = Parameters(8)*abs(sigmoid(t,Parameters(4),Parameters(6)) - sigmoid(t,Parameters(5),Parameters(6)));
            gammaf3  = Parameters(9)*sigmoid(t,Parameters(5),Parameters(6));
            gammaf   = gammaf1+gammaf2+gammaf3;
            dxdt    = zeros(3,1);                                    % x stands for "state" not x-axis, t is on x-axis
            dxdt(1) = -betag(i)*betaf/N*x(1)*x(2);                   % DE-1: S  (beta has now a new gain "betag")
            dxdt(2) = betag(i)*betaf/N*x(1)*x(2)-gammaf*x(2);        % DE-2: I  (beta has now a new gain "betag")
            dxdt(3) = gammaf*x(2);                                   % DE-3: R
            dxdt(4) = x(2);                                          % DE-4: C (you can eliminate this state by using "cumsum(I)' but for different step sizes, it might not be the best option)
            dSIR    = dxdt;
        end
        SIRout1 = SIRout0(:,2);                                      % Pick a state variable (here, we picked I)
    end
end