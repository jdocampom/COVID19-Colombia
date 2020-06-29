function f=COVID19(Name)
close all
 run(Name) % Read data from file 
Cases=Y;
L=length(Cases);
   
    figure('color','white')
    hold
    VP=[];
    L=length(Cases);
    Y=Cases;
    x=1:L;    
    plot(x,Cases,'.')   
    xx=1:300;
a=[ 10 10  .5 1]; % Set fitting parameters starting values
a=nlinfit(x,Cases,'corona2',a); % perform fit
Vpp=round(a(3)*atanh((((1+a(1)^2)^.5)-1)/a(1))+a(2));% Calculate inflection point
FF=corona2(a,xx); % Calculate fitted extrapolation
plot(xx,FF)
xlabel('Days Since Jan-22')
ylabel('Total No. of Infections')
title(Name)
grid
k=find(Cases<=.2*max(Cases)); % Select 80% of data
kk=max(k);
    U=[0 diff(Cases(kk:length(Y)))];% calculate daily increase (0 for arry length matching)
    E=max(U./Cases(kk:length(Y)));
UU=mean(U./Cases(kk:length(Y))); % Calculate mean daily tncreases
UUU=std(U./Cases(kk:length(Y)));% Caculate standard deviation of daily increases
CC=Cases;
Cases=CC;
Cases1=Cases; 
Cases2=Cases;
K=8; % Chose 8 extrapolated data points to be lifted up
for i=length(Y):length(Y)+K
       OLD1=(UU); % Setting lower limit 
       OLD2=(UU+3*UUU);% Setting upper limit to modefying eight point od data
      
    Cases1(i)=FF(i)*(1+OLD1); % Lifting up 8 points
    Cases2(i)=FF(i)*(1+OLD2);
  
    
   end

L1=length(Cases1);
   X=1:L1;
a1=a;
a2=a;
a1=nlinfit(X,Cases1,'corona2',a1); % Performing 
FF1=corona2(a1,xx);
plot(X(L1-8:L1),Cases1(L1-8:L1),'.k')% Plotting lifted up points
plot(xx,FF1,'--k')% Plotting new fit for lower limit
%______________________________
a2=nlinfit(X,Cases2,'corona2',a2);
FF2=corona2(a2,xx);
plot(X(L1-8:L1),Cases2(L1-8:L1),'.k')
plot(xx,FF2,'--k') % Plotting upper limit extrapolation


 U1=diff(FF1(Vpp:length(xx))); % Calculating daily changes to find saturation
 UU1=FF1(Vpp:length(xx)-1);
 ZZ1=U1./UU1;
 
 k1=find(ZZ1 <=0.001);% Finding Lower saturation
 if length(k1)==0
     k1=300;
 end
 %___________________________________
  U2=diff(FF2(Vpp:length(xx)));
 UU2=FF2(Vpp:length(xx)-1);
 ZZ2=U2./UU2;
 
 k2=find(ZZ2 <=0.001);% Finding Upper saturation
 if length(k2)==0
     k2=300;
 end
 

Name

Day1 = min(k1)+Vpp; % Date of Lower saturatio
Day2 = min(k2)+Vpp; % Date of Upper saturatio
Level1=max(FF1); % Magnitude of Lower saturation
Level2=max(FF2); % Magnitude of Upper saturation
Inflection_Point=Vpp
'After January-22-2020'
Upper_Limit=Level2
Lower_Limit=Level1
Upper_Date=Day2
'After January-22-2020'
Lower_Date=Day1
'After January-22-2020'
Summary=[Vpp Level2 Level1 Day2 Day1]% Arranging results
plot(Day1,corona2(a1,Day1),'or')% Plotting dates
plot(Day2,corona2(a2,Day2),'dr')
plot(Vpp,corona2(a,Vpp),'sr')
legend('Atual Data','Fit Using Actual Data Only','Upper Projected Points','Fit Using Upper Points', ...
    'Lower Projected Points','Fit Using Lower Points,','Saturartion Date 1', ...
    'Saturation Date 2', 'Inflection Point','location','southeast')
end
