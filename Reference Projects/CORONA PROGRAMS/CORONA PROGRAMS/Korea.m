function f=Koria
close all
Cases=[zeros(1,24 28 29 30 31 58 111 209 436 602 833 977 1261 1766 2337 3150 ...
    3736 4335 5186 5621 6284 6593 7041 7313 7478 7513 7755 7869 7979 ...
    8086 8162 8236 8320 8413 8565 8652 8799 8897 8961 9037 9137 9241 ...
    9332 9478 9583 9661 9786 9887 9976 10062 10156 10237 10284 10331 ...
    10384 10423 10450 10480 10512 10537 10564 10591 10613 10635 10653 ...
    10661 10674 10683 10694 10702 10708 10718 10728 10738 10752 10761 ...
    10765 10774 10780 10793 10801 10804 10806 10810 10822 10840 10874 ...
    10909 10939 10962 10991 11018 11037 11050 11065 11078 11110 11122 ...
    11142 11165 11190 11209 11225 11265 11344 11402 11441 11468 11503 ...
    11541 11590 11629 11668 11719 11776 11814 11852 11902 11947 12003 ...
    12051 12085 12121 12155 12198 12257];
    % Starting date 22 January End 18 June
Y=Cases;

% xx=1:length(Cases);
% U=Cases(length(Cases));
% W=max(xx);
% x=1:length(Cases);
% figure('color','white')
% plot(x,Cases,'k')
% xlabel('Days Since Jan-22')
% ylabel('Total No. of Infections')
% title('South Korea')
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
% if(LL>=.12)
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
% plot(W,U,'sr')
% plot(VP(1),corona2(b,VP(1)),'or')
% K=find(Y==0);
% Y(K)=[];
% Sumary=[max(corona2(AAA(1,:),250)),max(Y), floor(VP(1)),length(VP), M ,length(Y)]

% plot(A)
% xx=2:length(Cases);
% Cases=diff(Cases)
% UVx =[ 5
%     16
%     27
%     39
%     49
%     55
%     66
%     74
%     82
%     91
%     98
%    108
%    121
%    133
%    142
%    151
%    158
%    169
%    181
%    198];
% UV=[ 2.3697
%         2.654
%         3.128
%        3.6019
%        4.1706
%        4.8341
%        5.4028
%        6.1611
%        6.7299
%        7.2986
%        7.6777
%         8.436
%          8.91
%        9.2891
%        9.5735
%         9.763
%        9.9526
%        10.427
%        10.711
%        10.806
% ];
%    a1=[10 100 10];
%    b1=[10 100 10]'
%    a1=nlinfit(xx,Cases,'gauss',a1);
%    b1=nlinfit(UVx,UV,'gauss',b1);
%   
%    R=gauss(a1,xx)./gauss(b1,xx);
%    figure
%    length(xx)
%    length(Cases)
%    plot(xx,Cases)
%    hold
%    plot(xx,Gauss(a1,xx),'r')
%    figure
%    plot(UVx,UV)
%    hold
%    plot(xx,gauss(b1,xx),'r')
%    
%    a1
% 
% 
% 
% 
% 








