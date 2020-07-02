
%url=http://datosabiertos.esri.co/datasets/colombia-covid19-coronavirus-casos-diarios/data 
%% Read the excel file for covid cases (up to 22/05/20) 
data     = readtable('Colombia_COVID19.csv'); 
[row,col]= size(data); 
time     = [1:1:row]';                     %days until today 
n=length(time)*1e3;                                  %population 
infected_real   = table2array(data(:,3)); 
recovered_real  = table2array(data(:,5)); 
susceptible_real = n - infected_real-recovered_real; 

%% Computing beta and gamma 
b=n*[0 diff(susceptible_real)']./(-1.*(susceptible_real.*infected_real)); 
b(b==inf) =[]; 
b(b==-inf)=[]; 
b(isnan(b)) = []; 
b=reshape(b,[1,length(b)^2]); 
b=median(b) 
  
gamma=([0 diff(recovered_real)']./infected_real); 
gamma(gamma==inf) =[]; 
gamma(gamma==-inf)=[]; 
gamma(isnan(gamma)) = []; 
  
gamma=reshape(gamma,[1,length(gamma)^2]); 
gamma=median(gamma) 


%% Plotting model 
CI=[susceptible_real(1);1;0];    %FOR S,I,R 
time_span=[0 365*1]; 
[t,y] = ode45(@model,time_span,CI); 
  
plot(t,y(:,1),t,y(:,2),t,y(:,3)) 
legend(["Susceptible";"Infected";"Recovered"]) 
title('SIR modeling') 
xlabel('Days') 
ylabel('Population') 

function dydt = model(t,y) 
%y1 ->S(t) 
%y2 ->I(t) 
%y3 ->R(t) 
  
n=100000; 
beta=0.0766; 
gamma=0.0085; 
  
dydt = [-beta/n*y(1)*y(2);beta/n*y(1)*y(2)-gamma*y(2);gamma*y(2)]; 
end 
