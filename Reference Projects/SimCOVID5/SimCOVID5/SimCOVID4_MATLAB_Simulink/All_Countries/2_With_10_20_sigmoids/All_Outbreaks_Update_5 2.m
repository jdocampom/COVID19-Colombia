% SimCOVID Version 3.0 (MATLAB Program for Simulating COVID-19 Outbreak)
% Case: All World Countries
% Update: May 21, 2020
% Author: Ismael Abdulrahman
% Note 1: The data was taken from (updated daily): https://www.ecdc.europa.eu/en/geographical-distribution-2019-ncov-cases
% Note 2: You can control your curve-fitting by editing the initial values (Parameters0 vector).
% Note 3: The lower bounds for the parameters are equal or greater than zero except for the sigmoid parameter a. Changing the lower\upper bound limits can affect the solution and the time to solve the problem.
% Note 4: This program uses SIR model. It can be easily modified for new models
% Note 5: You only need to enter the country name @ line 12 and run
%==========================================================================
function Parameters = All_Outbreaks_Update_5
SearchCountry = 'italy';   % Enter the country name here (case insensitive). Use under-line dash "_" for a name with more than one part (ex: South_Africa, south_africa)
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
t00=findchangepts(Im,'Statistic','linear','MaxNumChanges',19);               % Find 19 points on x-axis that have major change in y (used for the initial values)
lt0=19-length(t00);
t0=[ones(1,lt0) t00]';
N0 = SelectedRowsFlipped.popData2018; N = N0(1);                             % Population
Im1=(smooth(Im))'; Im1(Im1==0)=1;
Im2=[0 Im1(1:end-1)];  Im3=Im2./Im1; Im4=normalize(Im3,'range',[0.1 1]);     % Find 9 initial values for the beta gains 2-10 that are close to confirmed cases (data)
Parameters000=Im4'; Parameters01=[Parameters000(t0(1)); Parameters000(t0) ];
Parameters0 = [Parameters01; t0; 0.1; 0.1]'; % Total initial values for 'lsqcurvefit'.
I0=Im(1); x0  = [N-I0; I0; 0; I0];                                           % Initial values for the states (S,I,R, C) used in 'ode45'
DateStart0 = SelectedRowsFlipped.dateRep; DateStart= DateStart0(1);          % First date with confirmed case(s)
Simulation_length = length(Im);                                              % Simulation time (day)
step_size = 0.1;                                                             % To give a smooth curve, choose a step size less than 1 (0.1 is preferred)
Simulated_t = 1:step_size:Simulation_length;
t_lb=[0; length(Cm)*1/19; length(Cm)*2/19; length(Cm)*3/19; length(Cm)*4/19;  length(Cm)*5/19; length(Cm)*6/19; length(Cm)*7/19; length(Cm)*8/19;  length(Cm)*9/19; length(Cm)*10/19; length(Cm)*11/19; length(Cm)*12/19;  length(Cm)*13/19; length(Cm)*14/19; length(Cm)*15/19; length(Cm)*16/19; length(Cm)*17/19; length(Cm)*18/19; ];              % lower band limits for the times used in the sigmoid function to be estimated
t_ub=[length(Cm)*1/19; length(Cm)*2/19; length(Cm)*3/19; length(Cm)*4/19;  length(Cm)*5/19; length(Cm)*6/19; length(Cm)*7/19; length(Cm)*8/19;  length(Cm)*9/19; length(Cm)*10/19; length(Cm)*11/19; length(Cm)*12/19;  length(Cm)*13/19; length(Cm)*14/19; length(Cm)*15/19; length(Cm)*16/19; length(Cm)*17/19; length(Cm)*18/19;  length(Cm)*19/19; ];  % upper band limits
%=============================Curve Fitting================================
[Parameters] = SIR_Model;
    function [Parameters] = SIR_Model
        options = optimoptions('lsqcurvefit','Algorithm','trust-region-reflective','Display', 'iter','FunctionTolerance',1e-16,'MaxFunctionEvaluations',10000*2,'MaxIterations',1000,'OptimalityTolerance',1e-16); % Increase iterations if the results are not suffiecent (algorithms: levenberg-marquardt, 'trust-region-reflective')
        %---------beta 1-20--------time 1-19  a   gamma
        lb = [ 0.01*(ones(1,20))'; t_lb;     -1;  0.01; ];      % minimum values for the parameters (a must). The sigmoid parameter (a) can be negative.
        ub = [ 10*(ones(1,20))';   t_ub;     10;  1  ];         % maximum values for the parameters (optional but preferred to be set)
        [Parameters] = lsqcurvefit(@SolveSIR1,Parameters0',tm',Cm', lb, ub, options);                      % Least-Square curve fitting
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
                plot(Simulated_t,Simulated_I, 'linewidth',2)
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
                plot(Simulated_t,Simulated_C, 'linewidth',2, 'Color',[0.9290 0.6940 0.1250])
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
            betaf1   = Parameters(1)*abs((sigmoid(Simulated_t,0,             Parameters(40)) - sigmoid(Simulated_t,Parameters(21),Parameters(40)))); % step function 1
            betaf2   = Parameters(2)*abs((sigmoid(Simulated_t,Parameters(21),Parameters(40)) - sigmoid(Simulated_t,Parameters(22),Parameters(40)))); % step function 2
            betaf3   = Parameters(3)*abs((sigmoid(Simulated_t,Parameters(22),Parameters(40)) - sigmoid(Simulated_t,Parameters(23),Parameters(40)))); % step function 3
            betaf4   = Parameters(4)*abs((sigmoid(Simulated_t,Parameters(23),Parameters(40)) - sigmoid(Simulated_t,Parameters(24),Parameters(40)))); % step function 4
            betaf5   = Parameters(5)*abs((sigmoid(Simulated_t,Parameters(24),Parameters(40)) - sigmoid(Simulated_t,Parameters(25),Parameters(40)))); % step function 2
            betaf6   = Parameters(6)*abs((sigmoid(Simulated_t,Parameters(25),Parameters(40)) - sigmoid(Simulated_t,Parameters(26),Parameters(40)))); % step function 3
            betaf7   = Parameters(7)*abs((sigmoid(Simulated_t,Parameters(26),Parameters(40)) - sigmoid(Simulated_t,Parameters(27),Parameters(40)))); % step function 4
            betaf8   = Parameters(8)*abs((sigmoid(Simulated_t,Parameters(27),Parameters(40)) - sigmoid(Simulated_t,Parameters(28),Parameters(40)))); % step function 2
            betaf9   = Parameters(9)*abs((sigmoid(Simulated_t,Parameters(28),Parameters(40)) - sigmoid(Simulated_t,Parameters(29),Parameters(40)))); % step function 3
            betaf10  = Parameters(10)*abs((sigmoid(Simulated_t,Parameters(29),Parameters(40)) - sigmoid(Simulated_t,Parameters(30),Parameters(40)))); % step function 2
            betaf11  = Parameters(11)*abs((sigmoid(Simulated_t,Parameters(30),Parameters(40)) - sigmoid(Simulated_t,Parameters(31),Parameters(40)))); % step function 3
            betaf12  = Parameters(12)*abs((sigmoid(Simulated_t,Parameters(31),Parameters(40)) - sigmoid(Simulated_t,Parameters(32),Parameters(40)))); % step function 4
            betaf13  = Parameters(13)*abs((sigmoid(Simulated_t,Parameters(32),Parameters(40)) - sigmoid(Simulated_t,Parameters(33),Parameters(40)))); % step function 2
            betaf14  = Parameters(14)*abs((sigmoid(Simulated_t,Parameters(33),Parameters(40)) - sigmoid(Simulated_t,Parameters(34),Parameters(40)))); % step function 3
            betaf15  = Parameters(15)*abs((sigmoid(Simulated_t,Parameters(34),Parameters(40)) - sigmoid(Simulated_t,Parameters(35),Parameters(40)))); % step function 2
            betaf16  = Parameters(16)*abs((sigmoid(Simulated_t,Parameters(35),Parameters(40)) - sigmoid(Simulated_t,Parameters(36),Parameters(40)))); % step function 3
            betaf17  = Parameters(17)*abs((sigmoid(Simulated_t,Parameters(36),Parameters(40)) - sigmoid(Simulated_t,Parameters(37),Parameters(40)))); % step function 4
            betaf18  = Parameters(18)*abs((sigmoid(Simulated_t,Parameters(37),Parameters(40)) - sigmoid(Simulated_t,Parameters(38),Parameters(40)))); % step function 2
            betaf19  = Parameters(19)*abs((sigmoid(Simulated_t,Parameters(38),Parameters(40)) - sigmoid(Simulated_t,Parameters(39),Parameters(40)))); % step function 3
            betaf20  = Parameters(20)*(sigmoid(Simulated_t,Parameters(39),Parameters(40)));
            betaf   = betaf1+betaf2+betaf3+betaf4+betaf5+ betaf6+betaf7+betaf8+betaf9+betaf10 +betaf11+betaf12+betaf13+betaf14+betaf15+ betaf16+betaf17+betaf18+betaf19+betaf20 ;
            gammaf=Parameters(41);
            
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
            betaf1   = Parameters(1)*abs((sigmoid(Simulated_t,0,             Parameters(40)) - sigmoid(Simulated_t,Parameters(21),Parameters(40)))); % step function 1
            betaf2   = Parameters(2)*abs((sigmoid(Simulated_t,Parameters(21),Parameters(40)) - sigmoid(Simulated_t,Parameters(22),Parameters(40)))); % step function 2
            betaf3   = Parameters(3)*abs((sigmoid(Simulated_t,Parameters(22),Parameters(40)) - sigmoid(Simulated_t,Parameters(23),Parameters(40)))); % step function 3
            betaf4   = Parameters(4)*abs((sigmoid(Simulated_t,Parameters(23),Parameters(40)) - sigmoid(Simulated_t,Parameters(24),Parameters(40)))); % step function 4
            betaf5   = Parameters(5)*abs((sigmoid(Simulated_t,Parameters(24),Parameters(40)) - sigmoid(Simulated_t,Parameters(25),Parameters(40)))); % step function 2
            betaf6   = Parameters(6)*abs((sigmoid(Simulated_t,Parameters(25),Parameters(40)) - sigmoid(Simulated_t,Parameters(26),Parameters(40)))); % step function 3
            betaf7   = Parameters(7)*abs((sigmoid(Simulated_t,Parameters(26),Parameters(40)) - sigmoid(Simulated_t,Parameters(27),Parameters(40)))); % step function 4
            betaf8   = Parameters(8)*abs((sigmoid(Simulated_t,Parameters(27),Parameters(40)) - sigmoid(Simulated_t,Parameters(28),Parameters(40)))); % step function 2
            betaf9   = Parameters(9)*abs((sigmoid(Simulated_t,Parameters(28),Parameters(40)) - sigmoid(Simulated_t,Parameters(29),Parameters(40)))); % step function 3
            betaf10  = Parameters(10)*abs((sigmoid(Simulated_t,Parameters(29),Parameters(40)) - sigmoid(Simulated_t,Parameters(30),Parameters(40)))); % step function 2
            betaf11  = Parameters(11)*abs((sigmoid(Simulated_t,Parameters(30),Parameters(40)) - sigmoid(Simulated_t,Parameters(31),Parameters(40)))); % step function 3
            betaf12  = Parameters(12)*abs((sigmoid(Simulated_t,Parameters(31),Parameters(40)) - sigmoid(Simulated_t,Parameters(32),Parameters(40)))); % step function 4
            betaf13  = Parameters(13)*abs((sigmoid(Simulated_t,Parameters(32),Parameters(40)) - sigmoid(Simulated_t,Parameters(33),Parameters(40)))); % step function 2
            betaf14  = Parameters(14)*abs((sigmoid(Simulated_t,Parameters(33),Parameters(40)) - sigmoid(Simulated_t,Parameters(34),Parameters(40)))); % step function 3
            betaf15  = Parameters(15)*abs((sigmoid(Simulated_t,Parameters(34),Parameters(40)) - sigmoid(Simulated_t,Parameters(35),Parameters(40)))); % step function 2
            betaf16  = Parameters(16)*abs((sigmoid(Simulated_t,Parameters(35),Parameters(40)) - sigmoid(Simulated_t,Parameters(36),Parameters(40)))); % step function 3
            betaf17  = Parameters(17)*abs((sigmoid(Simulated_t,Parameters(36),Parameters(40)) - sigmoid(Simulated_t,Parameters(37),Parameters(40)))); % step function 4
            betaf18  = Parameters(18)*abs((sigmoid(Simulated_t,Parameters(37),Parameters(40)) - sigmoid(Simulated_t,Parameters(38),Parameters(40)))); % step function 2
            betaf19  = Parameters(19)*abs((sigmoid(Simulated_t,Parameters(38),Parameters(40)) - sigmoid(Simulated_t,Parameters(39),Parameters(40)))); % step function 3
            betaf20  = Parameters(20)*(sigmoid(Simulated_t,Parameters(39),Parameters(40)));
            betaf   = betaf1+betaf2+betaf3+betaf4+betaf5+ betaf6+betaf7+betaf8+betaf9+betaf10 +betaf11+betaf12+betaf13+betaf14+betaf15+ betaf16+betaf17+betaf18+betaf19+betaf20 ;
            gammaf=Parameters(41); 
            
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
            betaf1   = Parameters(1)*abs((sigmoid(Simulated_t,0,             Parameters(40)) - sigmoid(Simulated_t,Parameters(21),Parameters(40)))); % step function 1
            betaf2   = Parameters(2)*abs((sigmoid(Simulated_t,Parameters(21),Parameters(40)) - sigmoid(Simulated_t,Parameters(22),Parameters(40)))); % step function 2
            betaf3   = Parameters(3)*abs((sigmoid(Simulated_t,Parameters(22),Parameters(40)) - sigmoid(Simulated_t,Parameters(23),Parameters(40)))); % step function 3
            betaf4   = Parameters(4)*abs((sigmoid(Simulated_t,Parameters(23),Parameters(40)) - sigmoid(Simulated_t,Parameters(24),Parameters(40)))); % step function 4
            betaf5   = Parameters(5)*abs((sigmoid(Simulated_t,Parameters(24),Parameters(40)) - sigmoid(Simulated_t,Parameters(25),Parameters(40)))); % step function 2
            betaf6   = Parameters(6)*abs((sigmoid(Simulated_t,Parameters(25),Parameters(40)) - sigmoid(Simulated_t,Parameters(26),Parameters(40)))); % step function 3
            betaf7   = Parameters(7)*abs((sigmoid(Simulated_t,Parameters(26),Parameters(40)) - sigmoid(Simulated_t,Parameters(27),Parameters(40)))); % step function 4
            betaf8   = Parameters(8)*abs((sigmoid(Simulated_t,Parameters(27),Parameters(40)) - sigmoid(Simulated_t,Parameters(28),Parameters(40)))); % step function 2
            betaf9   = Parameters(9)*abs((sigmoid(Simulated_t,Parameters(28),Parameters(40)) - sigmoid(Simulated_t,Parameters(29),Parameters(40)))); % step function 3
            betaf10  = Parameters(10)*abs((sigmoid(Simulated_t,Parameters(29),Parameters(40)) - sigmoid(Simulated_t,Parameters(30),Parameters(40)))); % step function 2
            betaf11  = Parameters(11)*abs((sigmoid(Simulated_t,Parameters(30),Parameters(40)) - sigmoid(Simulated_t,Parameters(31),Parameters(40)))); % step function 3
            betaf12  = Parameters(12)*abs((sigmoid(Simulated_t,Parameters(31),Parameters(40)) - sigmoid(Simulated_t,Parameters(32),Parameters(40)))); % step function 4
            betaf13  = Parameters(13)*abs((sigmoid(Simulated_t,Parameters(32),Parameters(40)) - sigmoid(Simulated_t,Parameters(33),Parameters(40)))); % step function 2
            betaf14  = Parameters(14)*abs((sigmoid(Simulated_t,Parameters(33),Parameters(40)) - sigmoid(Simulated_t,Parameters(34),Parameters(40)))); % step function 3
            betaf15  = Parameters(15)*abs((sigmoid(Simulated_t,Parameters(34),Parameters(40)) - sigmoid(Simulated_t,Parameters(35),Parameters(40)))); % step function 2
            betaf16  = Parameters(16)*abs((sigmoid(Simulated_t,Parameters(35),Parameters(40)) - sigmoid(Simulated_t,Parameters(36),Parameters(40)))); % step function 3
            betaf17  = Parameters(17)*abs((sigmoid(Simulated_t,Parameters(36),Parameters(40)) - sigmoid(Simulated_t,Parameters(37),Parameters(40)))); % step function 4
            betaf18  = Parameters(18)*abs((sigmoid(Simulated_t,Parameters(37),Parameters(40)) - sigmoid(Simulated_t,Parameters(38),Parameters(40)))); % step function 2
            betaf19  = Parameters(19)*abs((sigmoid(Simulated_t,Parameters(38),Parameters(40)) - sigmoid(Simulated_t,Parameters(39),Parameters(40)))); % step function 3
            betaf20  = Parameters(20)*(sigmoid(Simulated_t,Parameters(39),Parameters(40)));
            betaf   = betaf1+betaf2+betaf3+betaf4+betaf5+ betaf6+betaf7+betaf8+betaf9+betaf10 +betaf11+betaf12+betaf13+betaf14+betaf15+ betaf16+betaf17+betaf18+betaf19+betaf20 ;
            gammaf=Parameters(41);

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
