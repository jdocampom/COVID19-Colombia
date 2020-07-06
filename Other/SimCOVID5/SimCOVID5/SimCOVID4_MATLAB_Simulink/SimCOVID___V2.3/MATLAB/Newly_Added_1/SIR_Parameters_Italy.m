% SimCOVID Version 2.0, MATLAB Part, April 25, 2020
% Author: Ismael Abdulrahman
% The data was taken from (updated daily): https://www.ecdc.europa.eu/en/geographical-distribution-2019-ncov-cases
% Case: Italy Outbreak
% Note 1: You can control your curve-fitting by changing the values of st1, st2, and a. Choose critical points for st1 and st2 such as inflection points of the curve (Cm or Im). If not, choose st1 and st2 at one-third and two-third time of simulation (try it lastly). You can also add or remove some of the beta function branches. Choose "a" between 0 and 1 as initial values, if possible. Low values of "a" give a smoother curve.
% Note 2: The selected initial values are so so important. You can help the optimizer tool to find the solution much much faster if you choose these values wisely. 
% Note 3: The lower bound for the parameters here must give positive beta and gamma. Changing the upper bound limits can affect the solution and the time to solve the problem.
function [Parameters] = SIR_Parameters_Italy
I0  = 3;                            % Inital value of I (number of infectious befor simulation start time
R0  = 0;                            % Initial value of R "Recovered"
N   = 60*10^6;                      % Population
S0  = N-I0;                         % Initial value of "Susceptible"
st1 = 7; st2 = 22;                  % Step times for the sigmoid function (chosen empirically by observing Cm plot)
x0  = [S0; I0; R0 ];                % Initial values for the states (S,I,R) used in 'ode45'
[Parameters] = SIR_Model;
    function [Parameters] = SIR_Model
        I1  = xlsread('COVID-19-geographic-disbtribution-worldwide-2020-04-25','Italy','E:E');
        Im  = flip(I1);  Im=Im(54:end);            %  Neglecting zero elements and collecting total sum of infectious before simulation start time (I0=sum(Im(1:53)))
        Cm  = cumsum(Im);                          %  Cumulative sum of I
        tm  = 1:length(Cm);                        %  Time (day)
        %=============================Curve Fitting========================
        Parameters0 = [0.4; 0.32; 0.19; 1.44; 0.5];    % Initial values for 'lsqcurvefit'
        options = optimoptions('lsqcurvefit','Algorithm','trust-region-reflective','Display', 'iter','FunctionTolerance',1e-15,'MaxFunctionEvaluations',100000*2,'MaxIterations',10000,'OptimalityTolerance',1e-15); % levenberg-marquardt
        lb = [ 0; 0; 0; 0; 0];                       % minimum values for the parameters (a must)
        ub = 10*[ 1; 1; 1; 1; 1];                    % maximum values for the parameters (optional)
        [Parameters] = lsqcurvefit(@SolveSIR1,Parameters0,tm,Cm, lb, [], options);    % Least-Square curve fitting
        %=======================Printing SIR Parameters====================
        fprintf(1,'\tRate Constants:\n')
        for i = 1:length(Parameters)
            fprintf(1, '\t\tParameters(%d) = %8.5f\n', i, Parameters(i))
        end
        %==========================Plotting================================
        [~,Simulated_SIR] = SolveSIR2(Parameters, tm);
        Simulated_I  = Simulated_SIR(:,2);
        Simulated_S  = Simulated_SIR(:,1);
        Simulated_R  = Simulated_SIR(:,3);
        Simulated_C  = cumsum(Simulated_I);
        %===========================Plotting N=============================
        figure('Color',[1 1 1],'units','normalized')
        plot(tm,Simulated_S+Simulated_I+Simulated_R, 'linewidth',2)
        grid on;set(gca,'fontsize',16);
        xlabel('Time (day)','FontSize',16,'FontWeight','bold');
        xlim([0 length(tm)])
        ylabel('Total Population, N=S+I+R (person)','FontSize',16,'FontWeight','bold');
        dateaxis('x', 6, '22-Feb-2020')
        %===========================Plotting S=============================
        figure('Color',[1 1 1],'units','normalized')
        plot(tm,Simulated_S, 'linewidth',2)
        grid on;set(gca,'fontsize',16);
        xlabel('Time (day)','FontSize',16,'FontWeight','bold');
        xlim([0 length(tm)])
        ylabel('Susceptible (person)','FontSize',16,'FontWeight','bold');
        dateaxis('x', 6, '22-Feb-2020')
        %============================Plotting R============================
        figure('Color',[1 1 1],'units','normalized')
        plot(tm,Simulated_R, 'linewidth',2)
        grid on;set(gca,'fontsize',16);
        xlabel('Time (day)','FontSize',16,'FontWeight','bold');
        xlim([0 length(tm)])
        ylabel('Recovered (person)','FontSize',16,'FontWeight','bold');
        dateaxis('x', 6, '22-Feb-2020')
        %===========================Plotting C=============================
        figure('Color',[1 1 1],'units','normalized')
        plot(Simulated_C, 'linewidth',2)
        grid on;set(gca,'fontsize',16);
        xlabel('Time (day)','FontSize',16,'FontWeight','bold');
        xlim([0 length(tm)])
        ylabel('Cumulative Infectious (person)','FontSize',16,'FontWeight','bold');
        hold on
        plot(Cm,'o','LineWidth',1,'LineStyle',':','MarkerSize',6, 'Color','[0.8500, 0.3250, 0.0980]')
        legend({'simulated cumulative infectious','reported cases'}, 'FontSize',12);
        dateaxis('x', 6, '22-Feb-2020')
        %============================Plotting I============================
        figure('Color',[1 1 1],'units','normalized')
        plot(tm,Simulated_I, 'linewidth',2)
        grid on;set(gca,'fontsize',16);
        xlabel('Time (day)','FontSize',16,'FontWeight','bold');
        xlim([0 length(tm)])
        ylabel('Daily Infectious (person)','FontSize',16,'FontWeight','bold');
        hold on
        plot(Im,'p','LineWidth',1,'LineStyle',':','MarkerSize',10, 'Color','[0.8500, 0.3250, 0.0980]')
        legend({'simulated daily infectious','reported cases'}, 'FontSize',12);
        dateaxis('x', 6, '22-Feb-2020')
    end
%============================= Sigmoid Function============================
    function s = sigmoid(t,c,a)
        s = 1./(1 + exp(-a.*(t-c)));
    end
%============Solving SIR with One Output (Needed for lsqcurvefit)==========
    function SIRout1 = SolveSIR1(Parameters,t)
%         options = odeset('RelTol',1e-3,'AbsTol',1e-3);
        [~,SIRout0] = ode45(@DEs,t,x0);
        function dSIR=DEs(t,x)
            betaf1  = Parameters(1)*sigmoid(t,st1,-Parameters(4));  % See the Simulink file for this function
            betaf2  = Parameters(2)*abs(sigmoid(t,st1,Parameters(4)) - sigmoid(t,st2,Parameters(4)));
            betaf3  = Parameters(3)*sigmoid(t,st2,Parameters(4));
            betaf   = betaf1+betaf2+betaf3;
            dxdt    = zeros(3,1);                                   % x stands for "state" not x-axis, t is on x-axis
            dxdt(1) = -betaf/N*x(1)*x(2);                           % DE-1: S
            dxdt(2) = betaf/N*x(1)*x(2)-Parameters(5)*x(2);         % DE-2: I
            dxdt(3) = Parameters(5)*x(2);                           % DE-3: R
            dSIR    = dxdt;
        end
        SIRout1 = cumsum(SIRout0(:,2)) ;   % I is the second state of SIR. "cumsum" is for cumulative. We want to compare the cumulative infectious plots (reported and simulated) not the infectious variable
    end
%=================Solving SIR with All Outputs (for Plotting)==============
    function [tout, SIRout1] = SolveSIR2(Parameters,t)
%         options = odeset('RelTol',1e-3,'AbsTol',1e-3);
        [tout,SIRout0] = ode45(@DEs,t,x0);
        function dSIR=DEs(t,x)
            betaf1  = Parameters(1)*sigmoid(t,st1,-Parameters(4));  % See the Simulink file for this function
            betaf2  = Parameters(2)*abs(sigmoid(t,st1,Parameters(4)) - sigmoid(t,st2,Parameters(4)));
            betaf3  = Parameters(3)*sigmoid(t,st2,Parameters(4));
            betaf   = betaf1+betaf2+betaf3;
            dxdt    = zeros(3,1);                                   % x stands for "state" not x-axis, t is on x-axis
            dxdt(1) = -betaf/N*x(1)*x(2);                           % DE-1: S
            dxdt(2) = betaf/N*x(1)*x(2)-Parameters(5)*x(2);         % DE-2: I
            dxdt(3) = Parameters(5)*x(2);                           % DE-3: R
            dSIR    = dxdt;
        end
        SIRout1 = SIRout0;   
    end
end
