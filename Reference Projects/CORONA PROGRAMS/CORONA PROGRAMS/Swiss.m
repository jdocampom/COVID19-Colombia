
Cases=[zeros(1,24) 0 0 0 0 0 0 0 0 0 0 1 1 8 15 19 24 30 58 93 120 214 268 332 ...
    374 497 652 868 1139 1375 2217 2353 2742 3115 4222 5615 6863 7474 8795 ...
    9877 10897 11811 12928 14076 14829 15922 16606 17768 18827 19606 ...
    20505 21100 21657 22253 23280 24051 24551 25107 25415 25688 25936 ...
    26336 26732 27078 27404 27740 27944 28063 28268 28496 28677 28894 ...
    29061 29164 29264 29407 29586 29705 29817 29905 29981 30009 30060 ...
    30126 30207 30251 30305 30344 30380 30413 30463 30514 30572 30587 ...
    30597 30618 30658 30694 30707 30725 30736 30746 30761 30776 30796 ...
    30828 30845 30862 30871 30874 30893 30913 30936 30956 30965 30972 ...
    30988 31011 31044 31063 31094 31117 31131 31154 31187 31200];  % Starting date 22 Jan End 18 June 

Y=Cases;

% xx=1:length(Cases);
% U=Cases(length(Cases));
% W=max(xx);
% x=1:length(Cases);
% figure('color','white')
% plot(x,Cases,'k')
% xlabel('Days Since Jan-22')
% ylabel('Total No. of Infections')
% title('Switzerland')
% grid
% hold
% A=[];
% VP=[];
% AAA=[];
% for i=1:50
%     Cases=Cases(1:(length(x)-1*(i-1)));
%     xxx=1:length(Cases);
%    
% a=[ 10 10 .5 1];;
% a=nlinfit(xxx,Cases,'corona2',a);
% Vp=a(3)*atanh((((1+a(1)^2)^.5)-1)/a(1))+a(2);
% VP=[VP Vp];
% AA=a;
% A=[A a(1)];
% AAA=[AAA;a];
% L1=exp(min(A));
% L2=exp(max(A));
% LL=(L2-L1)/L2;
% if(LL>=.16)
%     M = length(Cases);
%     break
% else
% xx=1:150;
% plot(xx,corona2(a,xx),':k')
% end
% end
% b=AA(1,:)
% z=1:M;
% plot(z,Cases(1:M),'*k')
% plot(W,U,'sk')
% plot(VP(1),corona2(b,VP(1)),'or')
% K=find(Y==0);
% Y(K)=[];
% Sumary=[max(corona2(AAA(1,:),250)),max(Y), floor(VP(1)),length(VP), M ,length(Y)]

% Cases=diff(Cases);
% xx=1:length(Cases);
% UVx =[ 5
%     12
%     22
%     36
%     51
%     62
%     76
%     85
%     96
%    105
%    115
%    121
%    129
%    135
%    143
%    152
%    157
%    164
%    172
%    180];
% UV=[   0.65421
%       0.74766
%       0.93458
%         1.215
%        1.7757
%        2.5234
%         3.271
%        3.9252
%        4.5794
%        4.9533
%         5.514
%        5.8879
%        6.1682
%        6.5421
%        6.9159
%        7.2897
%        7.3832
%        7.6636
%        7.6636
%        7.9439];
%  a1=[10 100 10];
%    b1=[10 100 10];
%    a1=nlinfit(xx,Cases,'gauss',a1);
%    b1=nlinfit(UVx,UV,'gauss',b1);
%   
%    R=gauss(a1,xx)./gauss(b1,xx);
%    figure
%    length(xx)
%    length(Cases)
%    plot(xx,Cases)
%    hold
%    xxx=1:100;
%    plot(xxx,gauss(a1,xxx),'r')
%    figure
%    plot(UVx,UV)
%    hold
%    plot(xx,gauss(b1,xx),'r')
%    
%    
%    a1
%    
% 
% 
% 
% 




