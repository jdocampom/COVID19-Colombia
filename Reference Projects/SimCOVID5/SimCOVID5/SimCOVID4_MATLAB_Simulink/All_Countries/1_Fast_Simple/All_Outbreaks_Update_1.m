% SimCOVID Version 3.0 (MATLAB Program for COVID-19 Outbreak)
% Case: All World Countries
% Date: April 30, 2020
% Author: Ismael Abdulrahman
% Note 1: The data was taken from (updated daily): https://www.ecdc.europa.eu/en/geographical-distribution-2019-ncov-cases
% Note 2: You can control your curve-fitting by editing the initial values (Parameters0 vector).
% Note 3: The lower bounds for the parameters are equal or greater than zero. Changing the upper bound limits can affect the solution and the time to solve the problem.
% Note 4: This program uses SIR model. It can be easily modified for new models
% Note 5: You only need to enter the country name @ line 12 and run
%==========================================================================
function [Parameters] = All_Outbreaks_Update_1
SearchCountry = 'united_kingdom';   % Enter the country name here (case insensitive) and run. Use under-line dash "_" for a name with more than one part (ex: Saudi_Arabia, saudi_arabia)
COVID19_Data = readtable('COVID-19-geographic-disbtribution-worldwide-2020-05-15', 'ReadVariableNames', true);
RowIdx = find(strcmpi(SearchCountry, COVID19_Data.countriesAndTerritories)); % Row indices containing the input name
SelectedRows = COVID19_Data(RowIdx',:);                                      % Portion of the data with the given name
SelectedRowsFlipped = flipud(SelectedRows);                                  % Flip the rows (the original data is organized from new (top) to old (bottom))
SelectedRowsFlipped.CumulativeI = cumsum(SelectedRowsFlipped.cases);         % Add a column to the table to cumulatively sum I
toDelete = SelectedRowsFlipped.CumulativeI == 0;                             % Delete any zero cases before the first infectious is confirmed
SelectedRowsFlipped(toDelete,:) = [];                                        % The data after the deletion
Im = abs(SelectedRowsFlipped.cases);                                         % Measured (confirmed) infectious
Cm  = SelectedRowsFlipped.CumulativeI;                                       % Cumulative sum of Im
tm  = 1:length(Cm);                                                          % Time given with the data (day)
N0 = SelectedRowsFlipped.popData2018; N = N0(1);                             % Population
Parameters0 = [0.2; 0.3; 0.2; 0; length(Cm)/2; 0.1; 0.1];                    % Initial values for 'lsqcurvefit' (Any change is significantly affecting the solution). Ordered as: beta1, beta2, beta3, c1, c2, a, (sigmoid parameters), and gamma
I0=Im(1); x0  = [N-I0; I0; 0; I0];                                             % Initial values for the states (S,I,R) used in 'ode45'
DateStart0 = SelectedRowsFlipped.dateRep; DateStart= DateStart0(1);          % First date with confirmed cases
Simulation_length = length(Im);            % Simulation time (day)
step_size = 0.1;                    % To give a smooth curve, choose a step size less than 1 (0.1 is preferred)
Simulated_t = 1:step_size:Simulation_length;
%=============================Curve Fitting================================
[Parameters] = SIR_Model;
    function [Parameters] = SIR_Model
        options = optimoptions('lsqcurvefit','Algorithm','levenberg-marquardt','Display', 'iter','FunctionTolerance',1e-15,'MaxFunctionEvaluations',100000*2,'MaxIterations',2000,'OptimalityTolerance',1e-15); % Increase iterations if the results are not suffiecent (algorithms: levenberg-marquardt, 'trust-region-reflective')
        lb = [ 0; 0; 0; 0; 0; -1; 0];                                        % minimum values for the parameters (a must). The sigmoid parameter (a) can be negative
        ub = 1*[ 1; 1; 1; 1; 1; 1; 1];                                       % maximum values for the parameters (optional)
        [Parameters] = lsqcurvefit(@SolveSIR1,Parameters0,tm,Cm, lb, [], options);    % Least-Square curve fitting
        %=======================Printing SIR Parameters====================
        fprintf(1,'\tRate Constants:\n')
        for i = 1:length(Parameters)
            fprintf(1, '\t\tParameters(%d) = %8.5f\n', i, Parameters(i))
        end
        %==========================Plotting================================
        [~,Simulated_SIR] = SolveSIR2(Parameters, Simulated_t);
        Simulated_S  = Simulated_SIR(:,1);
        Simulated_I  = Simulated_SIR(:,2);
        Simulated_R  = Simulated_SIR(:,3);
        Simulated_C  = (Simulated_SIR(:,4));
        Simulated_N  = Simulated_S+Simulated_I+Simulated_R;
        %===========================Plotting N=============================
        figure('Color',[1 1 1],'units','normalized','outerposition',[0 0 1 1])
        plot(Simulated_t,Simulated_N, 'linewidth',2)
        grid on; grid minor; set(gca,'fontsize',16);
        xticks(0:Simulation_length/10:Simulation_length);
        xlabel('Time (day)','FontSize',16,'FontWeight','bold');
        xlim([0 Simulation_length])
        ylabel('N = S + I + R (person)','FontSize',16,'FontWeight','bold');
        dateaxis('x', 6, DateStart)
        %===========================Plotting S=============================
        figure('Color',[1 1 1],'units','normalized','outerposition',[0 0 1 1])
        plot(Simulated_t,Simulated_S, 'linewidth',2)
        grid on; grid minor; set(gca,'fontsize',16);
        xticks(0:Simulation_length/10:Simulation_length);
        xlabel('Time (day)','FontSize',16,'FontWeight','bold');
        xlim([0 Simulation_length])
        ylabel('Daily Susceptible (person)','FontSize',16,'FontWeight','bold');
        dateaxis('x', 6, DateStart)
        %===========================Plotting R=============================
        figure('Color',[1 1 1],'units','normalized','outerposition',[0 0 1 1])
        plot(Simulated_t,Simulated_R, 'linewidth',2)
        grid on; grid minor; set(gca,'fontsize',16);
        xticks(0:Simulation_length/10:Simulation_length);
        xlabel('Time (day)','FontSize',16,'FontWeight','bold');
        xlim([0 Simulation_length])
        ylabel('Daily Removed (person)','FontSize',16,'FontWeight','bold');
        dateaxis('x', 6, DateStart)
        %===========================Plotting C=============================
        figure('Color',[1 1 1],'units','normalized','outerposition',[0 0 1 1])
        plot(Simulated_t,Simulated_C, 'linewidth',2)
        grid on; grid minor; set(gca,'fontsize',16);
        xticks(0:Simulation_length/10:Simulation_length);
        xlabel('Time (day)','FontSize',16,'FontWeight','bold');
        xlim([0 Simulation_length])
        ylabel('Cumulative Infectious (person)','FontSize',16,'FontWeight','bold');
        hold on
        plot(Cm,'o','LineWidth',1,'LineStyle',':','MarkerSize',8, 'Color','[0.8500, 0.3250, 0.0980]')
        legend({'Cumulative infectious (simulated)','Daily infectious (reported cases'}, 'FontSize',12, 'location','northwest');
        dateaxis('x', 6, DateStart)
        %================================Plotting I========================
        figure('Color',[1 1 1],'units','normalized','outerposition',[0 0 1 1])
        plot(Simulated_t,Simulated_I, 'linewidth',2)
        grid on; grid minor; set(gca,'fontsize',16);
        xticks(0:Simulation_length/10:Simulation_length);
        xlabel('Time (day)','FontSize',16,'FontWeight','bold');
        xlim([0 Simulation_length])
        ylabel('Daily Infectious (person)','FontSize',16,'FontWeight','bold');
        hold on
        plot(Im,'p','LineWidth',1,'LineStyle',':','MarkerSize',8, 'Color','[0.8500, 0.3250, 0.0980]')
        legend({'Daily infectious (simulated)','Daily infectious (reported cases)'}, 'FontSize',12, 'location','northwest');
        dateaxis('x', 6, DateStart)
    end
%============================= Sigmoid Function============================
    function s = sigmoid(t,c,a)
        s = 1./(1 + exp(-a.*(t-c)));                                % This function is always positive regardless the signs of a and c
    end
%============Solving SIR with One Output (Needed for lsqcurvefit)==========
    function SIRout1 = SolveSIR1(Parameters,t)
        [~,SIRout0] = ode45(@DEs,t,x0);
        function dSIR=DEs(t,x)
            betaf1  = Parameters(1)*(sigmoid(t,0,           -Parameters(6)));
            betaf2  = Parameters(2)*abs((sigmoid(t,Parameters(4),Parameters(6)) - sigmoid(t,Parameters(5),Parameters(6))));
            betaf3  = Parameters(3)*(sigmoid(t,Parameters(5),Parameters(6)));
            betaf   = betaf1+betaf2+betaf3;
            dxdt    = zeros(4,1);                                   % x stands for "state" not x-axis, t is on x-axis
            dxdt(1) = -betaf/N*x(1)*x(2);                           % DE-1: S
            dxdt(2) = betaf/N*x(1)*x(2)-Parameters(7)*x(2);         % DE-2: I
            dxdt(3) = Parameters(7)*x(2);                           % DE-3: R
            dxdt(4) = x(2);                                         % DE-4: C (you can eliminate this state by using "cumsum' but for different step sizes, it might not be the best option)
            dSIR    = dxdt;
        end
        SIRout1 = cumsum(SIRout0(:,2)) ;
    end
%=================Solving SIR with All Outputs (for Plotting)==============
    function [tout, SIRout1] = SolveSIR2(Parameters,t)
        [tout,SIRout0] = ode45(@DEs,t,x0);
        function dSIR=DEs(t,x)
            betaf1  = Parameters(1)*(sigmoid(t,0,           -Parameters(6)));  % See the Simulink file for this function
            betaf2  = Parameters(2)*abs((sigmoid(t,Parameters(4),Parameters(6)) - sigmoid(t,Parameters(5),Parameters(6))));
            betaf3  = Parameters(3)*(sigmoid(t,Parameters(5),Parameters(6)));
            betaf   = betaf1+betaf2+betaf3;
            dxdt    = zeros(4,1);                                   % x stands for "state" not x-axis, t is on x-axis
            dxdt(1) = -betaf/N*x(1)*x(2);                           % DE-1: S
            dxdt(2) = betaf/N*x(1)*x(2)-Parameters(7)*x(2);         % DE-2: I
            dxdt(3) = Parameters(7)*x(2);                           % DE-3: R
            dxdt(4) = x(2);                                         % DE-4: C (you can eliminate this state by using "cumsum(I)' but for different step sizes, it might not be the best option)
            dSIR    = dxdt;
        end
        SIRout1 = SIRout0;                                          % Pick a state variable (here, we picked I)
    end
%==================================R0======================================
betaf1  = Parameters(1)*(sigmoid(Simulated_t,0,           -Parameters(6)));  % See the Simulink file for this function
betaf2  = Parameters(2)*abs((sigmoid(Simulated_t,Parameters(4),Parameters(6)) - sigmoid(Simulated_t,Parameters(5),Parameters(6))));
betaf3  = Parameters(3)*(sigmoid(Simulated_t,Parameters(5),Parameters(6)));
betaf   = betaf1+betaf2+betaf3;
R0=betaf'./Parameters(7);
figure('Color',[1 1 1],'units','normalized','outerposition',[0 0 1 1])
plot(Simulated_t,R0, 'linewidth',2)
grid on; grid minor; set(gca,'fontsize',16);
xticks(0:Simulation_length/10:Simulation_length);
xlabel('Time (day)','FontSize',16,'FontWeight','bold');
xlim([0 Simulation_length])
ylabel('Reproduction Number R_0 = \beta/\gamma','FontSize',16,'FontWeight','bold');
dateaxis('x', 6, DateStart)
end
%=================================END======================================
