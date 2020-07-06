% SimCOVID Version 4.0 (generalized)
% Case: All World Countries (tracking existing data and estimating future outbreaks)
% Date: May 25, 2020
% Author: Ismael Abdulrahman
% Note 1: The data was taken from (updated daily): https://www.ecdc.europa.eu/en/geographical-distribution-2019-ncov-cases
% Note 2: You can control your curve-fitting by editing (1) "MaxNumChanges" @line 34 and (2) the initial values (Parameters0 vector).
% Note 3: The lower bounds for the parameters are equal or greater than zero except for the sigmoid parameter a. Changing the lower\upper bound limits can significally affect the solution and the time to solve the problem.
% Note 4: This program uses the SIR model. It can be easily modified for new models
% Note 5: You only need to enter the country name @ line 18 and run the model
% Keynote 1: You can change the method of "findchangepts" @ line 34 from "rms" to "linear", "std", or "mean" and see the difference, especially for some cases (countries)
% Keynote 2: This is a generalized method for all outbreaks (current and future). For this reason, a smooth function @ line 37 is used to prevent some unexpected data. You can remove it for your case to see if that affects.
% Keynote 3: Change number of iterations (and other setting options) of the solver @line 65 for better results if needed
% Keynote 4: If the spread has not yet passed the peak (R0 > 1), you need to apply a pred-defined infection variable.
% Funding: The author has no funding source for this work (self-funded). This work has been continuously updated to support the community by the available personal sources.
%==========================================================================
function Parameters = Generalized_Method5_deaths
clc;
SearchCountry = 'united_states_of_america';   % Enter the country name here (case insensitive). Use under-line dash "_" for a name with more than one part (ex: South_Africa, south_africa)
COVID19_Data = readtable('COVID-19-geographic-disbtribution-worldwide-2020-05-25', 'ReadVariableNames', true);
RowIdx = find(strcmpi(SearchCountry, COVID19_Data.countriesAndTerritories)); % Row indices containing the input name
SelectedRows = COVID19_Data(RowIdx',:);                                      % Portion of the data extracted for the given country
SelectedRowsFlipped = flipud(SelectedRows);                                  % Flip the rows (the original data is organized from new (top) to old (bottom))
SelectedRowsFlipped.CumulativeI = cumsum(SelectedRowsFlipped.deaths);         % Add a column to the table to cumulatively sum I
toDelete = SelectedRowsFlipped.CumulativeI == 0;                             % Delete any zero cases before the first infectious is confirmed
SelectedRowsFlipped(toDelete,:) = [];                                        % The data after the deletion
Im = (abs(SelectedRowsFlipped.deaths))';                                      % Measured (confirmed) infectious
Cm = (SelectedRowsFlipped.CumulativeI)';                                     % Cumulative sum of Im
tm = 1:length(Cm);                                                           % Given time data (day)
%---------------------------Important Note---------------------------------
%---------------Number of sigmoid steps 'MaxNumChanges'--------------------
%------Keynote1: for tracking existing data, use high 'MaxNumChanges'------
%-----Keynote2: for estimating future outbreaks, vary 'MaxNumChanges'------
%-------------------Max. number of steps is length(Cm)---------------------
t0=findchangepts(Im,'Statistic','linear','MaxNumChanges',5);                % Find all points on x-axis (time) that have major change in y (it is also number of sigmoid functions). Set the value to any number such as 3, 5, 10,... or "length(Cm)". You can also enter t0 manually as a vector.
%--------------------End of number of sigmoid steps------------------------
N0 = SelectedRowsFlipped.popData2018; N = N0(1);                             % Population
Im1=(smooth(Im))'; Im1(Im1==0)=1;
Im2=[0 Im1(1:end-1)];  Im3=Im2./Im1;                                         % Find all initial values for the beta gains that are close to confirmed cases ratio (data) with a normalized range (use either "normalize" function or by dividing by max value
%-----------------------Normalize Im (two methods)-------------------------
% Im4=normalize(Im3,'range',[0.1 0.5]);                                      % Method 1 (change the range to get better results, if possible)
Im4=Im3./max(Im3);                                                           % Method 2
%---------------------------End of Normalize-------------------------------
Parameters000=Im4'; Parameters01= Parameters000(t0);
Parameters0 = [Parameters01; (t0(1:end-1))'; min(Parameters01); min(Parameters01)]';                    % Total initial values for 'lsqcurvefit'.
I0=Im(1); x0  = [N-I0; I0; 0; I0];                                           % Initial values for the states (S,I,R, C) used in 'ode45'
DateStart0 = SelectedRowsFlipped.dateRep; DateStart= DateStart0(1);          % First date with confirmed case(s)
%----------------------------Important Note--------------------------------
%--------------------------Simulation Length-------------------------------
%-----------Keynote1: use "length(Im)" for tracking existing data----------
%--------Keynore2: use > "length(Im)" for estimating future outbreaks------
%--------------Keynote3: for estimating future, use less t0----------------
Simulation_length = 1.5*length(Im);                                            % Simulation time (day)
%-------------------------End of Simulation Time---------------------------
step_size = 0.1;                                                             % To give a smooth curve, choose a step size less than 1 (0.1 is preferred)
Simulated_t = 1:step_size:Simulation_length;
t_b=((1:length(t0)-1)./length(t0)).*length(Cm);
t_lb=[0 (t_b(1:end-1)) ];                                                    % lower band limits for the times used in the sigmoid function to be estimated
t_ub=t_b;                                                                    % upper band limits
%-----------beta gains-----------time-----a-----gamma
lb = [ 0.01*ones(1,length(t0))   t_lb     0.25  0.01 ];                      % minimum values for the parameters (a must). The sigmoid parameter (a) can be negative.
ub = [ 10.0*ones(1,length(t0))   t_ub     0.75  1.00 ];                      % maximum values for the parameters (optional but preferred to be set)
%==========================Parameter Estimation============================
[Parameters] = SIR_Model;
    function [Parameters] = SIR_Model
        options = optimoptions('lsqcurvefit','Algorithm','trust-region-reflective','Display', 'iter','FunctionTolerance',1e-9,'MaxFunctionEvaluations',100000*2,'MaxIterations',300,'OptimalityTolerance',1e-9); % Increase iterations if the results are not suffiecent (algorithms: levenberg-marquardt, 'trust-region-reflective')
        %-------------------------Important Note---------------------------
        %--Compare either (1) Im (or gradient(Cm)) with simulated I (SIRout0(:,2)) or (2) Cm with simulated C (SIRout0(:,4))
        [Parameters] = lsqcurvefit(@SolveSIR1,Parameters0',tm',(gradient(Cm))', lb, ub, options);   % Least-Square curve fitting
        %=======================Printing SIR Parameters====================
        fprintf(1,'\tEstimated Parameters:\n')
        for s = 1:length(Parameters)
            fprintf(1, '\t\tParameters(%d) = %8.5f\n', s, Parameters(s))
        end
    end
%==========================Solving SIR Three Times=========================
Simulated_SIR=zeros(length(Simulated_t),3);
betag = [0.98 1.0 1.02];                                            % Beta change at time zero.
for k = 1:3
    [tout,Simulated_SIR(:,k)] = SolveSIR2(Parameters, Simulated_t); % Solve SIR with three values of beta (standard, plus minus a percentage)
end
%==========================Plotting all I's================================
t=tout;
I1 = Simulated_SIR(:,1); I2 = Simulated_SIR(:,2); I3 = Simulated_SIR(:,3);
tt = vertcat(t,flipud(t));
II = vertcat(I1,flipud(I3));
figure('Color',[1 1 1],'units','normalized','outerposition',[0 0 1 1])
bar(Im)
hold on
plot(t,I1, 'Color','[1 1 1]'); hold on
plot(t,I3, 'Color','[1 1 1]')
patch(tt,II,[0.75 0.75 0.75],'FaceAlpha',0.45, 'EdgeColor' , [1 1 1]) % trancperency (0 is for complete transp and 1 for fully opaque)
hold on
plot(t,I2,'b','LineWidth',3,'Color', '[0, 0.5, 0]')
% Use the color codes from this site: http://math.loyola.edu/~loberbro/matlab/html/colorsInMatlab.html
grid on;set(gca,'fontsize',16); grid minor;
xticks(0:Simulation_length/10:Simulation_length);
xlabel('Time (day)','FontSize',16,'FontWeight','bold');
xlim([0 Simulation_length])
ylabel('Daily Death Cases (person)','FontSize',16,'FontWeight','bold');
dateaxis('x', 6, DateStart)
%===========================Printing SIR Parameters========================
fprintf(1,'\tRate Constants:\n')
for i = 1:length(Parameters)
    fprintf(1, '\t\tParameters(%d) = %8.5f\n', i, Parameters(i))
end
%============================= Sigmoid Function============================
    function s = sigmoid(t,c,a)
        s = 1./(1 + exp(-a.*(t-c)));                                        % This function is always positive regardless the signs of a and c
    end
%============Solving SIR with One Output (Needed for lsqcurvefit)==========
    function SIRout1 = SolveSIR1(Parameters,Simulated_t)
        [~,SIRout0] = ode45(@DEs,Simulated_t,x0);
        function dSIR=DEs(Simulated_t,x)
            %--"sigmoid" is a smooth step-function proposed by the author for the beta (infection) function
            %--The time starts at zero, betaf1 is for this time
            %--Each step function has a starting point but infinte end, so we need to subtract two step-functions to obtain a step between two points (two times). The last time step does not have end time (betaf3).
            %--For speed reason, the for-loops are replaced by matrices (vectorization method, see betaf2)
            %--Sum of the sigmoid functions (betaf1,2,3)is the infection function
            betaf1   =     Parameters(1).*             abs((sigmoid(Simulated_t,0                                      ,Parameters(2*length(t0)) ))   -  (sigmoid(Simulated_t,Parameters(length(t0)+1)               ,Parameters(2*length(t0)) ))) ;
            betaf2   = sum(Parameters(2:length(t0)-1).*abs((sigmoid(Simulated_t,Parameters(length(t0)+1:2*length(t0)-2),Parameters(2*length(t0)) ))   -  (sigmoid(Simulated_t,Parameters(length(t0)+2:2*length(t0)-1),Parameters(2*length(t0)) ))));
            betaf3   =     Parameters(length(t0)    ).*   ((sigmoid(Simulated_t,Parameters(2*length(t0)-1)             ,Parameters(2*length(t0)) ))) ;
            betaf    = betaf1+betaf2+betaf3;
            gammaf=Parameters(2*length(t0)+1);
            %---------------------------DEs--------------------------------
            dxdt    = zeros(4,1);
            dxdt(1) = -betaf/N*x(1)*x(2);                            % DE-1: S
            dxdt(2) = betaf/N*x(1)*x(2)-gammaf*x(2);                 % DE-2: I
            dxdt(3) = gammaf*x(2);                                   % DE-3: R
            dxdt(4) = x(2);                                          % DE-4: C (you can eliminate this state by using "cumsum' or "cumtrapz" but for different step sizes, it might not be the best option)
            dSIR    = dxdt;
        end
        SIRout1 = (SIRout0(:,2)) ;                                   % Pick a state variable (here, I is selected to compare it with Im and minimize the difference)
    end
%===================Solving SIR Three Times for Plotting===================
    function [tout, SIRout1] = SolveSIR2(Parameters,Simulated_t)
        [tout,SIRout0] = ode45(@DEs,Simulated_t,x0);
        function dSIR=DEs(Simulated_t,x)
            betaf1   =     Parameters(1).*             abs((sigmoid(Simulated_t,0                                      ,Parameters(2*length(t0)) ))   -  (sigmoid(Simulated_t,Parameters(length(t0)+1)               ,Parameters(2*length(t0)) ))) ;
            betaf2   = sum(Parameters(2:length(t0)-1).*abs((sigmoid(Simulated_t,Parameters(length(t0)+1:2*length(t0)-2),Parameters(2*length(t0)) ))   -  (sigmoid(Simulated_t,Parameters(length(t0)+2:2*length(t0)-1),Parameters(2*length(t0)) ))));
            betaf3   =     Parameters(length(t0)    ).*   ((sigmoid(Simulated_t,Parameters(2*length(t0)-1)             ,Parameters(2*length(t0)) ))) ;
            betaf    = betaf1+betaf2+betaf3;
            gammaf=Parameters(2*length(t0)+1);
            %---------------------------DEs--------------------------------
            dxdt    = zeros(4,1);
            dxdt(1) = -betag(k)*betaf/N*x(1)*x(2);                            % DE-1: S
            dxdt(2) = betag(k)*betaf/N*x(1)*x(2)-gammaf*x(2);                 % DE-2: I
            dxdt(3) = gammaf*x(2);                                            % DE-3: R
            dxdt(4) = x(2);                                                   % DE-4: C (you can eliminate this state by using "cumsum' or "cumtrapz" but for different step sizes, it might not be the best option)
            dSIR    = dxdt;
        end
        SIRout1 = SIRout0(:,2);
    end
end
%=================================END======================================
