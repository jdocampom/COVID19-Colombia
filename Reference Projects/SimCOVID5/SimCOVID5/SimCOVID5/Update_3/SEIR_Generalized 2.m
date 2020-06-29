% SimCOVID Version 5.4 (Generalized SEIR Model - infection and recovery rates are variables)
% Case: All World Countries (tracking existing data and estimating future outbreaks)
% Date: June 03, 2020
% Author: Ismael Abdulrahman
% Note 1: The data was taken from (updated daily): https://www.ecdc.europa.eu/en/geographical-distribution-2019-ncov-cases
% Note 2: You can control your curve-fitting by editing (1) 'MaxNumChanges' and 'Statistic' @line 35 and (2) the initial values (Parameters0 vector) and their limits.
% Note 3: The lower bounds for the parameters are equal or greater than zero except for the sigmoid parameter "a". Changing the lower\upper bound limits can highly affect the solution and the time needed to solve the problem.
% Note 4: This program uses the SEIR model. It can be easily modified for new models
% Note 5: You only need to enter the country name @line-18 and run the model.
% Keynote 1: The number of sigmoid steps @line-35 plays a significant role especially in tracking data (high value gives smoother results and visa versa). For estimation, start with a low value of "MaxNumChanges" > 1, that is, 2, 3, 4, 5, and so on. You can increase it if you have enough data or the outbreak passed the peak. Typical values 3-5.
% Keynote 2: You can switch between the methods of "findchangepts" @ line 35 ("rms", "linear", "std", or "mean").
% Keynote 3: This is a generalized method for all outbreaks (current and future). For this reason, a "smooth" function @ line 38 is used to prevent some unexpected data. You can remove it for your case to see if that affects.
% Keynote 4: Change the number of iterations (and other setting options) of the solver @line-68 for better results if needed.
% Funding: The author has no funding source for this work (self-funded). This work has been continuously updated to support the community by the available personal resources.
%==========================================================================
function Parameters = SEIR_Generalized
clc;
SearchCountry = 'italy';   % Enter the country name here (case insensitive). Use under-line dash "_" for a name with more than one part (ex: South_Africa, south_africa)
COVID19_Data = readtable('COVID-19-geographic-disbtribution-worldwide-2020-06-01', 'ReadVariableNames', true);
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
t0=findchangepts(Im,'Statistic','mean','MaxNumChanges',5);                   % Find all points on x-axis (time) that have major change in y (it is also number of sigmoid functions). Set the value to any number such as 3, 5, 10,... or "length(Cm)". You can also enter t0 manually as a vector.
%--------------------End of number of sigmoid steps------------------------
N0 = SelectedRowsFlipped.popData2018; N = N0(1);                             % Population
Im1=(smooth(Im))'; Im1(Im1==0)=1;
Im2=[0 Im1(1:end-1)];  Im3=Im2./Im1;                                         % Find all initial values for the beta gains that are close to confirmed cases ratio (data) with a normalized range (use either "normalize" function or by dividing by max value
%-----------------------Normalize Im (two methods)-------------------------
% Im4=normalize(Im3,'range',[0.1 0.5]);                                      % Method 1 (change the range to get better results, if possible)
Im4=Im3./max(Im3);                                                           % Method 2
%---------------------------End of Normalize-------------------------------
Simulation_length = 1.75*length(Im);                                         % Simulation time (day). 
Parameters000=Im4'; Parameters01= Parameters000(t0);
incubation=normalize(1/rms(Parameters000),'range',[1/14 1/2]);
%--------------beta gains-----time-----------a--------------------expoent------------------------dc--lambda-gamma
Parameters0 = [Parameters01; (t0(1:end-1))'; min(Parameters000); mean(Parameters000)*length(Cm); 0; incubation;   min(Parameters000)]';       % Total initial values for 'lsqcurvefit'.
I0 = Im(1); S0 = N-I0;  E0 = I0; R0 = 0; x0  = [S0; E0; I0; R0];              % Initial values for the states (S,E,I,R) used in 'ode45'
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
%----------beta gains-----------time-----a--------expoent-----------dc------lambda-gamma
lb = [ 0.1*ones(1,length(t0))   t_lb     0.01     0.75*length(Cm)   0.0     1/14 0.010 ];              % minimum values for the parameters (a must). The sigmoid parameter (a) can be negative.
ub = [ 100*ones(1,length(t0))   t_ub     1.00     1.00*length(Cm)   0.1     1/2  1.00 ];               % maximum values for the parameters (optional but preferred to be set)
%==========================Parameter Estimation============================
[Parameters] = SEIR_Model;
    function [Parameters] = SEIR_Model
        options = optimoptions('lsqcurvefit','Algorithm','trust-region-reflective','Display', 'iter','FunctionTolerance',1e-9,'MaxFunctionEvaluations',100000*2,'MaxIterations',500,'OptimalityTolerance',1e-9); % Increase iterations if the results are not suffiecent (algorithms: levenberg-marquardt, 'trust-region-reflective')
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
betag = [0.975 1.0 1.025];                                                  % Beta change at times 1, 2, 3.
for k = 1:3
    [tout,Simulated_SEIR(:,k)] = SolveSEIR2(Parameters, Simulated_t);       % Solve SEIR with three values of beta (standard, plus minus a percentage)
end
%=============================Plotting C===================================
t=tout;
I1 = Simulated_SEIR(:,1); I2 = Simulated_SEIR(:,2); I3 = Simulated_SEIR(:,3);
tt = vertcat(t,flipud(t));
II = vertcat(I1,flipud(I3));
C1 = cumsum(Simulated_SEIR(:,1)); C2 = cumsum(Simulated_SEIR(:,2)); C3 = cumsum(Simulated_SEIR(:,3));
CC = vertcat(C1,flipud(C3));
figure('Color',[1 1 1],'units','normalized','outerposition',[0 0 1 1])      %----if you want to revesrse death and infectious stacked bars, use this: bar(tm, Im); hold on; bar(tm, Dm). If you want to plot only Im, use this: bar(Im)
bar(tm, Cm); hold on; bar(tm, cumsum(Dm))
hold on
plot(t,C1, 'Color','[1 1 1]'); hold on
plot(t,C3, 'Color','[1 1 1]')
patch(tt,CC,[0.75 0.75 0.75],'FaceAlpha',0.45, 'EdgeColor' , [1 1 1])       % trancperency (0 is for complete transp and 1 for fully opaque)
hold on
plot(t,C2,'b','LineWidth',3,'Color', '[0, 0.5, 0]')                         % Use the color codes from this site: http://math.loyola.edu/~loberbro/matlab/html/colorsInMatlab.html
grid on;set(gca,'fontsize',16); grid minor;
xticks(0:Simulation_length/10:Simulation_length);
xlabel('Time (day)','FontSize',16,'FontWeight','bold');
xlim([0 Simulation_length])
ylabel('Cumulative Infectious (person)','FontSize',16,'FontWeight','bold');
legend('cumulative infection cases', 'cumulative death cases','Location','southeast');
dateaxis('x', 6, DateStart)
%==========================Plotting all I's================================
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
legend('confirmed infection cases', 'confirmed death cases','Location','northeast');
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
            betaf3   =     Parameters(length(t0)  ).*     ((sigmoid(Simulated_t,Parameters(2*length(t0)-1)             ,Parameters(2*length(t0)) ))) ;
            betaf    = (betaf1+betaf2+betaf3)*exp(-Simulated_t/(Parameters(2*length(t0)+1))) + Parameters(2*length(t0)+2);   % Infection rate
            lambdaf  = Parameters(2*length(t0)+3);                          % Exposed rate
            gammaf   = Parameters(2*length(t0)+4);                          % Recovery rate
            %---------------------------DEs--------------------------------
            dxdt    = zeros(4,1);
            S = x(1); E = x(2); I = x(3); % R = x(4);
            dxdt(1) = -betaf/N*S*I;                                         % DE-1: S
            dxdt(2) = betaf/N*S*I-lambdaf*E;                                % DE-2: E
            dxdt(3) = lambdaf*E-gammaf*I;                                   % DE-3: I
            dxdt(4) = gammaf*I;                                             % DE-4: R
            dSEIR    = dxdt;
        end
        SEIRout1 = (SEIRout0(:,3)) ;                                        % Pick a state variable (here, I is selected to compare it with Im and minimize the difference)
    end
%===================Solving SEIR Three Times for Plotting==================
    function [tout, SEIRout1] = SolveSEIR2(Parameters,Simulated_t)
        [tout,SEIRout0] = ode45(@DEs,Simulated_t,x0);
        function dSEIR=DEs(Simulated_t,x)
            betaf1   =     Parameters(1).*             abs((sigmoid(Simulated_t,0                                      ,Parameters(2*length(t0)) ))   -  (sigmoid(Simulated_t,Parameters(length(t0)+1)               ,Parameters(2*length(t0)) ))) ;
            betaf2   = sum(Parameters(2:length(t0)-1).*abs((sigmoid(Simulated_t,Parameters(length(t0)+1:2*length(t0)-2),Parameters(2*length(t0)) ))   -  (sigmoid(Simulated_t,Parameters(length(t0)+2:2*length(t0)-1),Parameters(2*length(t0)) ))));
            betaf3   =     Parameters(length(t0)  ).*     ((sigmoid(Simulated_t,Parameters(2*length(t0)-1)             ,Parameters(2*length(t0)) ))) ;
            betaf    = (betaf1+betaf2+betaf3)*exp(-Simulated_t/(Parameters(2*length(t0)+1))) + Parameters(2*length(t0)+2);   % Infection rate
            lambdaf  = Parameters(2*length(t0)+3);                          % Exposed rate
            gammaf   = Parameters(2*length(t0)+4);                          % Recovery rate
            %---------------------------DEs--------------------------------
            dxdt    = zeros(4,1);
            S = x(1); E = x(2); I = x(3); % R = x(4);
            dxdt(1) = -betag(k)*betaf/N*S*I;                                % DE-1: S
            dxdt(2) = betag(k)*betaf/N*S*I-lambdaf*E;                       % DE-2: E
            dxdt(3) = lambdaf*E-gammaf*I;                                   % DE-3: I
            dxdt(4) = gammaf*I;                                             % DE-4: R
            dSEIR    = dxdt;
        end
        SEIRout1 = SEIRout0(:,3);
    end
% %===============================Plotting R0================================
% betaf1   =     Parameters(1).*             abs((sigmoid(Simulated_t,0                                      ,Parameters(2*length(t0)) ))   -  (sigmoid(Simulated_t,Parameters(length(t0)+1)               ,Parameters(2*length(t0)) ))) ;
% betaf2   = sum(((Parameters(2:length(t0)-1)).*abs((sigmoid(Simulated_t,Parameters(length(t0)+1:2*length(t0)-2),Parameters(2*length(t0)) ))   -  (sigmoid(Simulated_t,Parameters(length(t0)+2:2*length(t0)-1),Parameters(2*length(t0)) )))));
% betaf3   =     Parameters(length(t0)    ).*   ((sigmoid(Simulated_t,Parameters(2*length(t0)-1)             ,Parameters(2*length(t0)) ))) ;
% betaf    = (betaf1+betaf2+betaf3).*exp(-Simulated_t/(Parameters(2*length(t0)+1))) + Parameters(2*length(t0)+2);   % Infection rate
% lambdaf  = Parameters(2*length(t0)+3);                          % Exposed rate
% gammaf   = Parameters(2*length(t0)+4);                          % Recovery rate
% R0=betaf./gammaf;
% figure('Color',[1 1 1],'units','normalized','outerposition',[0 0 1 1])
% plot(Simulated_t,R0, 'linewidth',2);
% grid on; grid minor; set(gca,'fontsize',16);
% xticks(0:Simulation_length/10:Simulation_length);
% xlabel('Time (day)','FontSize',16,'FontWeight','bold');
% xlim([0 Simulation_length])
% ylabel('Reproduction Number R_0 = \beta/\gamma','FontSize',16,'FontWeight','bold');
% dateaxis('x', 6, DateStart)
end
%=================================END======================================
