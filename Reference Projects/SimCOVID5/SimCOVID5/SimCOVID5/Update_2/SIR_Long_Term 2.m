% SimCOVID Version 5.0 (generalized SEIR Model, Long-Term)
% Case: All World Countries (tracking existing data and estimating future outbreaks)
% Date: May 31, 2020
% Author: Ismael Abdulrahman
% Note 1: The data was taken from (updated daily): https://www.ecdc.europa.eu/en/geographical-distribution-2019-ncov-cases
% Note 2: You can control your curve-fitting by editing (1) 'MaxNumChanges' and 'Statistic' @line 35 and (2) the initial values (Parameters0 vector).
% Note 3: The lower bounds for the parameters are equal or greater than zero except for the sigmoid parameter a. Changing the lower\upper bound limits can significally affect the solution and the time to solve the problem.
% Note 4: This program uses the SEIR model. It can be easily modified for new models
% Note 5: You only need to enter the country name @ line 18 and run the model
% Keynote 1: You can change the method of "findchangepts" @ line 35 from "rms" to "linear", "std", or "mean" and see the difference, especially for some cases (countries)
% Keynote 2: This is a generalized method for all outbreaks (current and future). For this reason, a smooth function @ line 38 is used to prevent some unexpected data. You can remove it for your case to see if that affects.
% Keynote 3: Change number of iterations (and other setting options) of the solver @line 66 for better results if needed
% Keynote 4: If the spread has not yet passed the peak (R0 > 1), you need to apply a pred-defined infection variable.
% Funding: The author has no funding source for this work (self-funded). This work has been continuously updated to support the community by the available personal sources.
%==========================================================================
function Parameters = SIR_Long_Term
clc;
SearchCountry = 'italy';   % Enter the country name here (case insensitive). Use under-line dash "_" for a name with more than one part (ex: South_Africa, south_africa)
COVID19_Data = readtable('COVID-19-geographic-disbtribution-worldwide-2020-05-31', 'ReadVariableNames', true);
RowIdx = find(strcmpi(SearchCountry, COVID19_Data.countriesAndTerritories)); % Row indices containing the input name
SelectedRows = COVID19_Data(RowIdx',:);                                      % Portion of the data extracted for the given country
SelectedRowsFlipped = flipud(SelectedRows);                                  % Flip the rows (the original data is organized from new (top) to old (bottom))
SelectedRowsFlipped.CumulativeI = cumsum(SelectedRowsFlipped.cases);         % Add a column to the table to cumulatively sum I
toDelete = SelectedRowsFlipped.CumulativeI == 0;                             % Delete any zero cases before the first infectious is confirmed
SelectedRowsFlipped(toDelete,:) = [];                                        % The data after the deletion
Im = (abs(SelectedRowsFlipped.cases))';                                      % Measured (confirmed) infectious
Dm = (abs(SelectedRowsFlipped.deaths))';                                     % Measured (confirmed) deads
Cm = (SelectedRowsFlipped.CumulativeI)';                                     % Cumulative sum of Im
tm = 1:length(Cm);                                                           % Given time data (day)
%---------------------------Important Note---------------------------------
%---------------Number of sigmoid steps 'MaxNumChanges'--------------------
%------Keynote1: for tracking existing data, use high 'MaxNumChanges'------
%-----Keynote2: for estimating future outbreaks, vary 'MaxNumChanges'------
%-------------------Max. number of steps is length(Cm)---------------------
t0=findchangepts(Im,'Statistic','linear','MaxNumChanges',5);                 % Find all points on x-axis (time) that have major change in y (it is also number of sigmoid functions). Set the value to any number such as 3, 5, 10,... or "length(Cm)". You can also enter t0 manually as a vector.
%--------------------End of number of sigmoid steps------------------------
N0 = SelectedRowsFlipped.popData2018; N = N0(1);                             % Population
Im1=(smooth(Im))'; Im1(Im1==0)=1;
Im2=[0 Im1(1:end-1)];  Im3=Im2./Im1;                                         % Find all initial values for the beta gains that are close to confirmed cases ratio (data) with a normalized range (use either "normalize" function or by dividing by max value
%-----------------------Normalize Im (two methods)-------------------------
% Im4=normalize(Im3,'range',[0.1 0.5]);                                      % Method 1 (change the range to get better results, if possible)
Im4=Im3./max(Im3);                                                           % Method 2
%---------------------------End of Normalize-------------------------------
Simulation_length = 1.75*length(Im);                                         % Simulation time (day). Note: use this carefully as it is now part in the beta function.
Parameters000=Im4'; Parameters01= Parameters000(t0);
Parameters0 = [Parameters01; (t0(1:end-1))'; min(Parameters01); Simulation_length*1; 0; min(Parameters01)]';       % Total initial values for 'lsqcurvefit'.
I0 = Im(1); S0 = N-I0;  R0 = 0; x0  = [S0; I0; R0];              % Initial values for the states (S,E,I,R) used in 'ode45'
DateStart0 = SelectedRowsFlipped.dateRep; DateStart= DateStart0(1);          % First date with confirmed case(s)
%----------------------------Important Note--------------------------------
%--------------------------Simulation Length-------------------------------
%-----------Keynote1: use "length(Im)" for tracking existing data----------
%--------Keynore2: use > "length(Im)" for estimating future outbreaks------
%--------------Keynote3: for estimating future, use less t0----------------                                       
%-------------------------End of Simulation Time---------------------------
step_size = 1;                                                               
Simulated_t = 1:step_size:Simulation_length;
t_b=((1:length(t0)-1)./length(t0)).*length(Cm);
t_lb=[0 (t_b(1:end-1)) ];                                                    % lower band limits for the times used in the sigmoid function to be estimated
t_ub=t_b;                                                                    % upper band limits
%-----------beta gains-----------time-----a-----lambda-gamma
lb = [ 0.01*ones(1,length(t0))   t_lb     0.25   Simulation_length*1/4   0     0.01 ];               % minimum values for the parameters (a must). The sigmoid parameter (a) can be negative.
ub = [ 10.0*ones(1,length(t0))   t_ub     0.75   Simulation_length*4/4   0.1   10.0 ];               % maximum values for the parameters (optional but preferred to be set)
%==========================Parameter Estimation============================
[Parameters] = SEIR_Model;
    function [Parameters] = SEIR_Model
        options = optimoptions('lsqcurvefit','Algorithm','trust-region-reflective','Display', 'iter','FunctionTolerance',1e-9,'MaxFunctionEvaluations',100000*2,'MaxIterations',1000,'OptimalityTolerance',1e-9); % Increase iterations if the results are not suffiecent (algorithms: levenberg-marquardt, 'trust-region-reflective')
        %-------------------------Important Note---------------------------
        %--Compare either (1) Im (or gradient(Cm)) with simulated I (SEIRout0(:,2)) or (2) Cm with simulated C (SEIRout0(:,4))
        [Parameters] = lsqcurvefit(@SolveSEIR1,Parameters0',tm',(gradient(Cm))', lb, ub, options);   % Least-Square curve fitting
        %=======================Printing SEIR Parameters===================
        fprintf(1,'\tEstimated Parameters:\n')
        for s = 1:length(Parameters)
            fprintf(1, '\t\tParameters(%d) = %8.5f\n', s, Parameters(s))
        end
    end
%==========================Solving SEIR Three Times========================
Simulated_SEIR=zeros(length(Simulated_t),3);
betag = [0.975 1.0 1.025];                                                  % Beta change at time zero.
for k = 1:3
    [tout,Simulated_SEIR(:,k)] = SolveSEIR2(Parameters, Simulated_t);       % Solve SEIR with three values of beta (standard, plus minus a percentage)
end
%==========================Plotting all I's================================
t=tout;
I1 = Simulated_SEIR(:,1); I2 = Simulated_SEIR(:,2); I3 = Simulated_SEIR(:,3);
tt = vertcat(t,flipud(t));
II = vertcat(I1,flipud(I3));
figure('Color',[1 1 1],'units','normalized','outerposition',[0 0 1 1])
bar(tm, [Im; Dm],'stacked');                                                %----if you want to revesrse death and infectious stacked bars, use this: bar(tm, Im); hold on; bar(tm, Dm). If you want to plot only Im, use this: bar(Im)
hold on
plot(t,I1, 'Color','[1 1 1]'); hold on
plot(t,I3, 'Color','[1 1 1]')
patch(tt,II,[0.75 0.75 0.75],'FaceAlpha',0.45, 'EdgeColor' , [1 1 1])       % trancperency (0 is for complete transp and 1 for fully opaque)
hold on
plot(t,I2,'b','LineWidth',3,'Color', '[0, 0.5, 0]')                         % Use the color codes from this site: http://math.loyola.edu/~loberbro/matlab/html/colorsInMatlab.html
grid on;set(gca,'fontsize',16); grid minor;
xticks(0:Simulation_length/10:Simulation_length);
xlabel('Time (day)','FontSize',16,'FontWeight','bold');
xlim([0 Simulation_length])
ylabel('Daily Infectious (person)','FontSize',16,'FontWeight','bold');
legend('confirmed infection cases', 'confirmed death cases');
dateaxis('x', 6, DateStart)
%===========================Printing SEIR Parameters=======================
fprintf(1,'\tRate Constants:\n')
for i = 1:length(Parameters)
    fprintf(1, '\t\tParameters(%d) = %8.5f\n', i, Parameters(i))
end
%============================= Sigmoid Function============================
    function s = sigmoid(t,c,a)
        s = 1./(1 + exp(-a.*(t-c)));                                        % This function is always positive regardless the signs of a and c
    end
%============Solving SEIR with One Output (Needed for lsqcurvefit)=========
    function SEIRout1 = SolveSEIR1(Parameters,Simulated_t)
        [~,SEIRout0] = ode45(@DEs,Simulated_t,x0);
        function dSEIR=DEs(Simulated_t,x)
            %--"sigmoid" is a smooth step-function proposed by the author for the beta (infection) function
            %--The time starts at zero, betaf1 is for this time
            %--Each step function has a starting point but infinte end, so we need to subtract two step-functions to obtain a step between two points (two times). The last time step does not have end time (betaf3).
            %--For speed reason, the for-loops are replaced by matrices (vectorization method, see betaf2)
            %--Sum of the sigmoid functions (betaf1,2,3)is the infection function
            betaf1   =     Parameters(1).*             abs((sigmoid(Simulated_t,0                                      ,Parameters(2*length(t0)) ))   -  (sigmoid(Simulated_t,Parameters(length(t0)+1)               ,Parameters(2*length(t0)) ))) ;
            betaf2   = sum(Parameters(2:length(t0)-1).*abs((sigmoid(Simulated_t,Parameters(length(t0)+1:2*length(t0)-2),Parameters(2*length(t0)) ))   -  (sigmoid(Simulated_t,Parameters(length(t0)+2:2*length(t0)-1),Parameters(2*length(t0)) ))));
            betaf3   =     Parameters(2*length(t0)  ).*   ((sigmoid(Simulated_t,Parameters(2*length(t0)-1)             ,Parameters(2*length(t0)) ))) ;
            betaf    = (betaf1+betaf2+betaf3)*abs(1-Simulated_t/Parameters(2*length(t0)+1)) + Parameters(2*length(t0)+2);   % Infection rate
            gammaf   = Parameters(2*length(t0)+3);                          % Recovery rate
            %---------------------------DEs--------------------------------
            dxdt    = zeros(3,1);
            S = x(1); I = x(2); % R = x(3); 
            dxdt(1) = -betaf/N*S*I;                                         % DE-1: S
            dxdt(2) = betaf/N*S*I-gammaf*I;                                 % DE-2: I
            dxdt(3) = gammaf*I;                                             % DE-3: R
            dSEIR   = dxdt;
        end
        SEIRout1 = (SEIRout0(:,2)) ;                                        % Pick a state variable (here, I is selected to compare it with Im and minimize the difference)
    end
%===================Solving SEIR Three Times for Plotting==================
    function [tout, SEIRout1] = SolveSEIR2(Parameters,Simulated_t)
        [tout,SEIRout0] = ode45(@DEs,Simulated_t,x0);
        function dSEIR=DEs(Simulated_t,x)
            betaf1   =     Parameters(1).*             abs((sigmoid(Simulated_t,0                                      ,Parameters(2*length(t0)) ))   -  (sigmoid(Simulated_t,Parameters(length(t0)+1)               ,Parameters(2*length(t0)) ))) ;
            betaf2   = sum(Parameters(2:length(t0)-1).*abs((sigmoid(Simulated_t,Parameters(length(t0)+1:2*length(t0)-2),Parameters(2*length(t0)) ))   -  (sigmoid(Simulated_t,Parameters(length(t0)+2:2*length(t0)-1),Parameters(2*length(t0)) ))));
            betaf3   =     Parameters(2*length(t0)  ).*   ((sigmoid(Simulated_t,Parameters(2*length(t0)-1)             ,Parameters(2*length(t0)) ))) ;
            betaf    = (betaf1+betaf2+betaf3)*abs(1-Simulated_t/Parameters(2*length(t0)+1)) + Parameters(2*length(t0)+2);   % Infection rate
            gammaf   = Parameters(2*length(t0)+3);                          % Recovery rate
            %---------------------------DEs--------------------------------
            dxdt    = zeros(3,1);
            S = x(1); I = x(2); % R = x(3); 
            dxdt(1) = -betag(k)*betaf/N*S*I;                                % DE-1: S
            dxdt(2) = betag(k)*betaf/N*S*I-gammaf*I;                        % DE-2: I
            dxdt(3) = gammaf*I;                                             % DE-3: R
            dSEIR   = dxdt;
        end
        SEIRout1 = SEIRout0(:,2);
    end
end
%=================================END======================================
