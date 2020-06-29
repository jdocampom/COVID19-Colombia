% SimCOVID Version 4.0 (generalized)
% Case: All World Countries
% Date: May 21, 2020
% Author: Ismael Abdulrahman
% Note 1: The data was taken from (updated daily): https://www.ecdc.europa.eu/en/geographical-distribution-2019-ncov-cases
% Note 2: You can control your curve-fitting by editing the initial values (Parameters0 vector).
% Note 3: The lower bounds for the parameters are equal or greater than zero except for the sigmoid parameter a. Changing the lower\upper bound limits can affect the solution and the time to solve the problem.
% Note 4: This program uses the SIR model. It can be easily modified for new models
% Note 5: You only need to enter the country name @ line 17 and run
% Keynote 1: You can change the method of "findchangepts" @ line 28 from "rms" to "linear", "std", or "mean" and see the difference, especially for some cases (countries)
% Keynote 2: This is a generalized method for all outbreaks (current and future). For this reason, a smooth function @ line 30 is used to prevent some unexpected data. You can remove it for your case to see if that affects.
% Keynote 3: Change number of iterations (and other setting options) of the solver @line 48 for better results if needed 
% Funding: The author has no funding source for this work (self-funded). This work has been continuously updated to support the community by available personal means.
%==========================================================================
function Parameters = Generalized_Method
clc;
SearchCountry = 'iran';   % Enter the country name here (case insensitive). Use under-line dash "_" for a name with more than one part (ex: South_Africa, south_africa)
COVID19_Data = readtable('COVID-19-geographic-disbtribution-worldwide-2020-05-21', 'ReadVariableNames', true);
RowIdx = find(strcmpi(SearchCountry, COVID19_Data.countriesAndTerritories)); % Row indices containing the input name
SelectedRows = COVID19_Data(RowIdx',:);                                      % Portion of the data extracted for the given country
SelectedRowsFlipped = flipud(SelectedRows);                                  % Flip the rows (the original data is organized from new (top) to old (bottom))
SelectedRowsFlipped.CumulativeI = cumsum(SelectedRowsFlipped.cases);         % Add a column to the table to cumulatively sum I
toDelete = SelectedRowsFlipped.CumulativeI == 0;                             % Delete any zero cases before the first infectious is confirmed
SelectedRowsFlipped(toDelete,:) = [];                                        % The data after the deletion
Im = (abs(SelectedRowsFlipped.cases))';                                      % Measured (confirmed) infectious
Cm  = (SelectedRowsFlipped.CumulativeI)';                                    % Cumulative sum of Im
tm  = 1:length(Cm);                                                          % Time given with the data (day)
t00=findchangepts(Im,'Statistic','rms','MaxNumChanges',length(Cm));          % Find all points on x-axis that have major change in y (used for the initial values)
N0 = SelectedRowsFlipped.popData2018; N = N0(1);                             % Population
Im1=(smooth(Im))'; Im1(Im1==0)=1;
Im2=[0 Im1(1:end-1)];  Im3=Im2./Im1; Im4=normalize(Im3,'range',[0.1 1]);     % Find all initial values for the beta gains that are close to confirmed cases ratio (data) with a normalized range
Parameters000=Im4'; Parameters01= Parameters000(t00);
Parameters0 = [Parameters01; (t00(1:end-1))'; min(Parameters01); min(Parameters01)]';                    % Total initial values for 'lsqcurvefit'.
I0=Im(1); x0  = [N-I0; I0; 0; I0];                                           % Initial values for the states (S,I,R, C) used in 'ode45'
DateStart0 = SelectedRowsFlipped.dateRep; DateStart= DateStart0(1);          % First date with confirmed case(s)
Simulation_length = length(Im);                                              % Simulation time (day)
step_size = 0.1;                                                             % To give a smooth curve, choose a step size less than 1 (0.1 is preferred)
Simulated_t = 1:step_size:Simulation_length;
t_b=((1:length(t00)-1)./length(t00)).*length(Cm);
t_lb=[0 (t_b(1:end-1)) ];                                                    % lower band limits for the times used in the sigmoid function to be estimated
t_ub=t_b;                                                                    % upper band limits
%-----------beta gains------------time-----a-----gamma
lb = [ 0.01*ones(1,length(t00))   t_lb     -1    0.01 ];                     % minimum values for the parameters (a must). The sigmoid parameter (a) can be negative.
ub = [ 10.0*ones(1,length(t00))   t_ub     10    1    ];                     % maximum values for the parameters (optional but preferred to be set)
%=============================Curve Fitting================================
[Parameters] = SIR_Model;
    function [Parameters] = SIR_Model
        options = optimoptions('lsqcurvefit','Algorithm','trust-region-reflective','Display', 'iter','FunctionTolerance',1e-6,'MaxFunctionEvaluations',10000*2,'MaxIterations',200,'OptimalityTolerance',1e-6); % Increase iterations if the results are not suffiecent (algorithms: levenberg-marquardt, 'trust-region-reflective')
        [Parameters] = lsqcurvefit(@SolveSIR1,Parameters0',tm',Cm', lb, ub, options);   % Least-Square curve fitting
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
        ylabel('Susceptible (person)','FontSize',16,'FontWeight','bold');
        dateaxis('x', 6, DateStart)
        %===========================Plotting R=============================
        figure('Color',[1 1 1],'units','normalized','outerposition',[0 0 1 1])
        plot(Simulated_t,Simulated_R, 'linewidth',2)
        grid on; grid minor; set(gca,'fontsize',16);
        xticks(0:Simulation_length/10:Simulation_length);
        xlabel('Time (day)','FontSize',16,'FontWeight','bold');
        xlim([0 Simulation_length])
        ylabel('Removed (person)','FontSize',16,'FontWeight','bold');
        dateaxis('x', 6, DateStart)
        %==========================Plotting I (bar)=======================
        figure('Color',[1 1 1],'units','normalized','outerposition',[0 0 1 1])
        bar(Im,'FaceColor','flat','LineWidth',0.5)
        hold on
        plot(Simulated_t,Simulated_I, 'linewidth',4)
        grid on; grid minor; set(gca,'fontsize',16);
        xticks(0:Simulation_length/10:Simulation_length);
        xlabel('Time (day)','FontSize',16,'FontWeight','bold');
        xlim([0 Simulation_length])
        ylabel('Daily Infectious (person)','FontSize',16,'FontWeight','bold');
        legend({'Daily infectious (reported cases)','Daily infectious (simulated)'}, 'FontSize',12, 'location','northwest');
        dateaxis('x', 6, DateStart)
        %==========================Plotting C (bar)========================
        figure('Color',[1 1 1],'units','normalized','outerposition',[0 0 1 1])
        bar(Cm,'FaceColor','flat','LineWidth',0.5)
        hold on
        plot(Simulated_t,Simulated_C, 'linewidth',4)
        grid on; grid minor; set(gca,'fontsize',16);
        xticks(0:Simulation_length/10:Simulation_length);
        xlabel('Time (day)','FontSize',16,'FontWeight','bold');
        xlim([0 Simulation_length])
        ylabel('Cumulative Infectious (person)','FontSize',16,'FontWeight','bold');
        legend({'Cumulative infectious (reported cases)','Cumulative infectious (simulated)'}, 'FontSize',12, 'location','northwest');
        dateaxis('x', 6, DateStart)
        %================================Plotting I========================
        figure('Color',[1 1 1],'units','normalized','outerposition',[0 0 1 1])
        plot(Im,'o','LineWidth',1,'LineStyle',':','MarkerSize',8)
        hold on
        plot(Simulated_t,Simulated_I, 'linewidth',2)
        grid on; grid minor; set(gca,'fontsize',16);
        xticks(0:Simulation_length/10:Simulation_length);
        xlabel('Time (day)','FontSize',16,'FontWeight','bold');
        xlim([0 Simulation_length])
        ylabel('Daily Infectious (person)','FontSize',16,'FontWeight','bold');
        legend({'Daily infectious (reported cases)','Daily infectious (simulated)'}, 'FontSize',12, 'location','northwest');
        dateaxis('x', 6, DateStart)
        %===========================Plotting C=============================
        figure('Color',[1 1 1],'units','normalized','outerposition',[0 0 1 1])
        plot(Cm,'o','LineWidth',1,'LineStyle',':','MarkerSize',8)
        hold on
        plot(Simulated_t,Simulated_C, 'linewidth',2)
        grid on; grid minor; set(gca,'fontsize',16);
        xticks(0:Simulation_length/10:Simulation_length);
        xlabel('Time (day)','FontSize',16,'FontWeight','bold');
        xlim([0 Simulation_length])
        ylabel('Cumulative Infectious (person)','FontSize',16,'FontWeight','bold');
        legend({'Cumulative infectious (reported cases)','Cumulative infectious (simulated)'}, 'FontSize',12, 'location','northwest');
        dateaxis('x', 6, DateStart)
    end
%============================= Sigmoid Function============================
    function s = sigmoid(t,c,a)
        s = 1./(1 + exp(-a.*(t-c)));                                        % This function is always positive regardless the signs of a and c
    end
%============Solving SIR with One Output (Needed for lsqcurvefit)==========
    function SIRout1 = SolveSIR1(Parameters,Simulated_t)
        [~,SIRout0] = ode45(@DEs,Simulated_t,x0);
        function dSIR=DEs(Simulated_t,x)
            betaf1   =     Parameters(1).*              abs((sigmoid(Simulated_t,0                                        ,Parameters(2*length(t00)) ))   -  (sigmoid(Simulated_t,Parameters(length(t00)+1)                ,Parameters(2*length(t00)) ))) ;  % First function starts at time zero.
            betaf2   = sum(Parameters(2:length(t00)-1).*abs((sigmoid(Simulated_t,Parameters(length(t00)+1:2*length(t00)-2),Parameters(2*length(t00)) ))   -  (sigmoid(Simulated_t,Parameters(length(t00)+2:2*length(t00)-1),Parameters(2*length(t00)) ))) ); % The sigmoid functions between the first and last one
            betaf3   =     Parameters(length(t00)    ).*   ((sigmoid(Simulated_t,Parameters(2*length(t00)-1)              ,Parameters(2*length(t00)) ))) ;
            betaf    = betaf1+betaf2+betaf3;
            gammaf=Parameters(2*length(t00)+1);
            
            dxdt    = zeros(4,1);
            dxdt(1) = -betaf/N*x(1)*x(2);                            % DE-1: S
            dxdt(2) = betaf/N*x(1)*x(2)-gammaf*x(2);                 % DE-2: I
            dxdt(3) = gammaf*x(2);                                   % DE-3: R
            dxdt(4) = x(2);                                          % DE-4: C (you can eliminate this state by using "cumsum' or "cumtrapz" but for different step sizes, it might not be the best option)
            dSIR    = dxdt;
        end
        SIRout1 = cumsum(SIRout0(:,2)) ;                             % Pick a state variable (here, I is selected to compare it with Im and minimize the difference)
    end
%=================Solving SIR with All Outputs (for Plotting)==============
    function [tout, SIRout1] = SolveSIR2(Parameters,Simulated_t)
        [tout,SIRout0] = ode45(@DEs,Simulated_t,x0);
        function dSIR=DEs(Simulated_t,x)
            betaf1   =     Parameters(1).*              abs((sigmoid(Simulated_t,0                                        ,Parameters(2*length(t00)) ))   -  (sigmoid(Simulated_t,Parameters(length(t00)+1),Parameters(2*length(t00)) ))) ;
            betaf2   = sum(Parameters(2:length(t00)-1).*abs((sigmoid(Simulated_t,Parameters(length(t00)+1:2*length(t00)-2),Parameters(2*length(t00)) ))   -  (sigmoid(Simulated_t,Parameters(length(t00)+2:2*length(t00)-1),Parameters(2*length(t00)) ))) );
            betaf3   =     Parameters(length(t00)    ).*   ((sigmoid(Simulated_t,Parameters(2*length(t00)-1)                ,Parameters(2*length(t00)) ))) ;
            betaf    = betaf1+betaf2+betaf3;
            gammaf=Parameters(2*length(t00)+1);
            
            dxdt    = zeros(4,1);
            dxdt(1) = -betaf/N*x(1)*x(2);                            % DE-1: S
            dxdt(2) = betaf/N*x(1)*x(2)-gammaf*x(2);                 % DE-2: I
            dxdt(3) = gammaf*x(2);                                   % DE-3: R
            dxdt(4) = x(2);                                          % DE-4: C (you can eliminate this state by using "cumsum' or "cumtrapz" but for different step sizes, it might not be the best option)
            dSIR    = dxdt;
        end
        SIRout1 = SIRout0;
    end
%===============================Plotting R0================================
betaf1   =     Parameters(1).*              abs((sigmoid(Simulated_t,0                                        ,Parameters(2*length(t00)) ))   -  (sigmoid(Simulated_t,Parameters(length(t00)+1),Parameters(2*length(t00)) ))) ;
betaf2   = sum(Parameters(2:length(t00)-1).*abs((sigmoid(Simulated_t,Parameters(length(t00)+1:2*length(t00)-2),Parameters(2*length(t00)) ))   -  (sigmoid(Simulated_t,Parameters(length(t00)+2:2*length(t00)-1),Parameters(2*length(t00)) ))));
betaf3   =     Parameters(length(t00)    ).*   ((sigmoid(Simulated_t,Parameters(2*length(t00)-1)              ,Parameters(2*length(t00)) ))) ;
betaf    = betaf1+betaf2+betaf3;
gammaf=Parameters(2*length(t00)+1);
R0=betaf./gammaf;

figure('Color',[1 1 1],'units','normalized','outerposition',[0 0 1 1])
plot(Simulated_t,R0, 'linewidth',2);
grid on; grid minor; set(gca,'fontsize',16);
xticks(0:Simulation_length/10:Simulation_length);
xlabel('Time (day)','FontSize',16,'FontWeight','bold');
xlim([0 Simulation_length])
ylabel('Reproduction Number R_0 = \beta/\gamma','FontSize',16,'FontWeight','bold');
dateaxis('x', 6, DateStart)
end
%=================================END======================================
