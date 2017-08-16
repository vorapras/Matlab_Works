%% Code for FM457A and FM457B - Lecture 6

% clear workspace
clear all;
close all;
clc;


%% No Aggregate Risk
%draw random numbers and simulation
N=500; T=1000;
mui=ones(1,N)*1.15;
sigmai=.1+.55*rand(1,N);

%form returns, build portfolios by consecutively adding a security, and
%plot it
for i=1:N;
    R(:,i)=mui(i)+sigmai(i)*randn(T,1);
end;
for i=2:N;
    P=mean(R(:,1:i)')';
    portstd(i)=std(P);
end;
plot(portstd,'b'); hold on;


%% With Aggregate Risk
%draw random numbers and simulation
N=500; T=1000;
mui=ones(1,N)*1.15;
sigmai=.1+.25*rand(1,N);
betaM=ones(1,N); Y=.16*randn(T,1);

for i=1:N;
    RA(:,i)=mui(i)+sigmai(i)*randn(T,1)+betaM(i)*Y;
end;
for i=2:500;
    P=mean(RA(:,1:i)')';
    portstdA(i)=std(P);
end;
plot(portstdA, 'r'); hold off;

%% One Factor Model
%draw random numbers and simulate factor (market) return
T=1000; mu=1.1; sigma=.16; x=randn(T,1); 
rM=mu+sigma*x;
rf=1.02;

%simulate cross-section
N=10000;
B=-.5+2.5*rand(1,N);
sigmai=.05+.55*rand(1,N);
mui=rf+B*(mu-rf);
R=crosssection(mu,sigma,rf,B,sigmai,T);

%create portfolios
A=[2 5 10 20 50 100 500];
for j=1:7;
  for i=1:22;
    x(i)=(.98+.01*i);
    s=0;
    in=zeros(1,N);
    for k=1:N;
      if s<=A(j) & abs(mean(R(:,k))-x(i))<.005;
        s=s+1;
        in(k)=1;
      end;
    end;
    P=mean(R(:,in==1)')';
    portmean(i,j)=mean(P); %portfolio mean return
    portstd(i,j)=std(P);   %portfolio risk
  end;
end;

%plot mean-stdev frontier
plot(portstd(:,1),portmean(:,1),'c');
hold on;
plot(portstd(:,2),portmean(:,2),'r--');
plot(portstd(:,3),portmean(:,3),'b');
plot(portstd(:,4),portmean(:,4),'m--');
plot(portstd(:,5),portmean(:,5),'g');
plot(portstd(:,6),portmean(:,6),'k--');
plot(portstd(:,7),portmean(:,7),'y');
plot(0,rf,'rd');
legend('2','5','10','20','50','100','500','rf');
