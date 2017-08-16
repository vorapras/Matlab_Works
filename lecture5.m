%% Code for FM457A and FM457B - Lecture 5

% clear workspace
clear all;
close all;
clc;


%% Simulating GARCH

%%% Discuss function simGARCH %%%

%generate simulated GARCH data
w=.05; a=.1; b=.8; T=1000;
out=simGARCH(w,a,b,T);
subplot(3,1,1); plot(out(:,1)); title('r');
subplot(3,1,2); plot(out(:,1).^2); title('r^{2}');
subplot(3,1,3); plot(out(:,2)); title('\sigma^{2}');

%% Estimating GARCH Manually

%generate simulated GARCH data
w=.05; a=.1; b=.8; T=10000;
out=simGARCH(w,a,b,T);

%estimate on simulated data
n=100; clear X;
for i=1:n;
    X(:,i)=out(n-i+1:T-i,1).^2;
end;
Y=out(n+1:T,1).^2;
regcoef=regress(Y,[ones(T-n,1) X]);
aest=regcoef(2);
best=regcoef(3)/regcoef(2);
west=regcoef(1)/((1-best^n)/(1-best));
disp([w a b;west aest best]);

%check how accurate different equations hold
subplot;
plot(a*b.^[0:8],regcoef(2:10),'.');
hold on;
plot(regcoef(2:10),regcoef(2:10),'k');
xlabel('Coefficients implied by true parameters');
ylabel('Actual coefficients');
%all points should line up along the 45 degree line


%% Estimate GARCH Using MATLAB Routines

% For new econometric toolbox
% ToEstMdl = garch(1,1);
% EstMdl= estimate(ToEstMdl,out(:,1));

[w_est, a_est, b_est] = ugarch(out(:,1) , 1 , 1);

disp([w a b]);
disp([w_est a_est b_est]);

%% Compare Simulated with Actual Data

%generate simulated Microsoft returns and plot along actual returns
data_msft;
rmsft=msft(:,4);
sigma=std(rmsft); mu=mean(rmsft)-.5*sigma^2; T=length(rmsft);
timeline=round(msft(1,2)/10000)+(round(msft(T,2)/10000)-round(msft(1,2)/10000))*([1:T]-1)/(T-1);
nlow=sum(rmsft<-.1);
nhigh=sum(rmsft>.1);
p1=nlow/T; p2=nhigh/T; 
J1=-.15; J2=.15;
sigmaS=.0025; rho=.96;
rsim=simsecSV(mu,sigma,J1,J2,p1,p2,rho,sigmaS,T);
subplot(2,1,1); plot(timeline,rmsft(1:T)); 
title('MSFT'); axis([timeline(1) timeline(T) -.2 .2]);
subplot(2,1,2); plot(timeline,rsim(1:T)); 
title('Simulated'); axis([timeline(1) timeline(T) -.2 .2]);

%%
%estimate GARCH(1,1) manually on actual data
n=100; clear X;
for i=1:n;
    X(:,i)=rmsft(n-i+1:T-i,1).^2;
end;
Y=rmsft(n+1:T,1).^2;
regcoef=regress(Y,[ones(T-n,1) X]);
a=regcoef(2);
b=regcoef(3)/regcoef(2);
w=regcoef(1)/((1-b^n)/(1-b));
disp([w a b]);

%estimate GARCH(1,1) manually on simulated data
n=100; clear X;
for i=1:n;
    X(:,i)=rsim(n-i+1:T-i,1).^2;
end;
Y=rsim(n+1:T,1).^2;
regcoefS=regress(Y,[ones(T-n,1) X]);
asim=regcoefS(2);
bsim=regcoefS(3)/regcoefS(2);
wsim=regcoefS(1)/((1-bsim^n)/(1-bsim));
disp([wsim asim bsim]);

%test how well the unused equations hold
subplot;
plot(a*b.^[0:8],regcoef(2:10),'sb','MarkerSize',10);
hold on;
plot(asim*bsim.^[0:8],regcoefS(2:10),'ro','MarkerSize',10);
plot(regcoef(2:10),regcoef(2:10),'k');
xlabel('Coefficients implied by \alpha, \beta');
ylabel('Actual coefficients');
legend('MSFT','Simulated');

%% Estimate Volatility Time Series from Actual Microsoft Returns
vmsft=zeros(T,1);
n=100;
for t=n+1:T;
    k=0; s=0;
    for i=1:n;
        k=k+b^(i-1);
        s=s+(rmsft(t-i)^2)*b^(i-1); %make use of equations on slide 50
    end;
    vmsft(t)=sqrt(w*k+a*s);
end;
vmsft(1:n)=mean(vmsft(n+1:T));
subplot(2,1,1); plot(timeline,rmsft(1:T));
title('MSFT Return'); axis([timeline(1) timeline(T) -.2 .2]);
subplot(2,1,2); plot(timeline,vmsft(1:T));
title('MSFT Volatility'); axis([timeline(1) timeline(T) 0 .06]);

%% Value at Risk - Definition

T0=100000; sigma=.023;
x=randn(T0,1)*sigma;
rvar5=sqrt(2*sigma^2)*erfinv(2*.05-1);
hist(x,50); axis([-.15 .15 0 8000]); hold on;
plot(rvar5*ones(6,1), 0:8000/5:8000,'r','LineWidth',3);
disp([rvar5 sum(x<rvar5)/T0]);

%% Value at Risk - Back Testing

%In order to run this part you need to load MSFT data, estimate MSFT
%GARCH(1,1) parameters, and estimate MSFT volatility (has been done above)
var5msft=sqrt(2*vmsft.^2)*erfinv(2*.05-1);
subplot(2,1,1); plot(timeline,vmsft);
title('MSFT Volatility'); axis([timeline(1) timeline(T) 0 .06]);
subplot(2,1,2); 
plot(timeline,rmsft,'b'); hold on;
plot(timeline,var5msft,'r','LineWidth',2);
title('MSFT VAR 5%'); axis([timeline(1) timeline(T) -.15 .05]);

%plot volatility, returns and VaR to see the violations
disp('Percent Violations GARCH VAR 5%');
disp(sum(rmsft<var5msft)/T);
var5msftconstv=sqrt(2*std(rmsft).^2)*erfinv(2*.05-1);
disp('Percent Violations Constant Vol VAR 5%');
disp(sum(rmsft<var5msftconstv)/T);
%Note that both overestimate MSFT's number of extreme returns and therefore
%MSFT's risk


        