% SimCOVID Version 2.0, MATLAB Part, April 25, 2020
% Author: Ismael Abdulrahman
% The data was taken from (updated daily): https://www.ecdc.europa.eu/en/geographical-distribution-2019-ncov-cases
% Case: USA Outbreak
function [Parameters] = SIR_Parameters_USA
N   = 328*10^6;                       % Population
st1 = 22; st2 = 26;                  % Step times for the sigmoid function (chosen empirically by observing Cm plot)
x0  = [N/1; 16; 0; ];                % Initial values for the states (S,I,R) used in 'ode45'
[Parameters] = SIR_Model;
    function [Parameters] = SIR_Model
        I1  = xlsread('COVID-19-geographic-disbtribution-worldwide-2020-04-25','USA','E:E');
        Im  = flip(I1);  Im=Im(53:end);            %  Neglecting zero elements considering 14 starting-cases (I0=14)
        Cm  = cumsum(Im);                          %  Cumulative sum of I
        tm  = 1:length(Cm);                        %  Time (day)
        %=============================Curve Fitting========================
        Parameters0 = [0.4; 0.3; 0.1; 0.25; 0.2];    % Initial values for 'lsqcurvefit'
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
        Simulated_C  = cumsum(Simulated_I);
        %===========================Plotting C=============================
        figure('Color',[1 1 1],'units','normalized','outerposition',[0 0 1 1])
        plot(Simulated_C, 'linewidth',3)
        grid on;set(gca,'fontsize',16);
        xlabel('Time (day)','FontSize',16,'FontWeight','bold');
        xlim([0 length(tm)])
        ylabel('Cumulative Infectious (person)','FontSize',16,'FontWeight','bold');
        hold on
        plot(Cm,'*','LineWidth',1,'LineStyle',':','MarkerSize',10, 'Color','[0.8500, 0.3250, 0.0980]')
        legend({'simulated cumulative infectious','reported cases'}, 'FontSize',12);
        dateaxis('x', 6, '21-Feb-2020')
        %================================Plotting I========================
        figure('Color',[1 1 1],'units','normalized','outerposition',[0 0 1 1])
        plot(tm,Simulated_I, 'linewidth',3)
        grid on;set(gca,'fontsize',16);
        xlabel('Time (day)','FontSize',16,'FontWeight','bold');
        xlim([0 length(tm)])
        ylabel('Daily Infectious (person)','FontSize',16,'FontWeight','bold');
        hold on
        plot(Im,'*','LineWidth',1,'LineStyle',':','MarkerSize',10, 'Color','[0.8500, 0.3250, 0.0980]')
        legend({'simulated daily infectious','reported cases'}, 'FontSize',12);
        dateaxis('x', 6, '21-Feb-2020')
    end
%============================= Sigmoid Function============================
    function s = sigmoid(t,c,a)
        s = 1./(1 + exp(-a.*(t-c)));
    end
%============Solving SIR with One Output (Needed for lsqcurvefit)==========
    function SIRout1 = SolveSIR1(Parameters,t)
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
