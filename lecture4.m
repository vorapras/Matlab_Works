%% Code for FM457A and FM457B - Lecture 4

% clear workspace
clear all;
close all;
clc;


%% Simulating Jumps

%load data
data_msft;

%%% Discuss Function simsec %%%

%calibrate jump probabilities
subplot(2,1,1); 
hist(msft(:,4),[-.2:.01:.2]);
axis([-.2 .2 0 800]); xlabel('Actual');
sigma=std(msft(:,4)); mu=mean(msft(:,4))-.5*sigma^2; T=length(msft); 
nlow=sum(msft(:,4)<-.1);
nhigh=sum(msft(:,4)>.1);
[T a]=size(msft);
p1=nlow/T; p2=nhigh/T; 
J1=-.15; J2=.15;

%%% Discuss Function simsecJ %%%

%simulate Microsoft returns
r=simsecJ(mu,sigma,J1,J2,p1,p2,T);
disp([mean(msft(:,4)) mean(r)]);
disp([std(msft(:,4)) std(r)]);
disp([skewness(msft(:,4)) skewness(r)]);
disp([kurtosis(msft(:,4)) kurtosis(r)]);
subplot(2,1,2); hist(r,[-.2:.01:.2]); 
axis([-.2 .2 0 800]); xlabel('Simulated');
%Note that kurtosis of simulated now matches actual


%% Volatility

%moving average
T1=floor(T/7); ma=zeros(T1,2);
for i=1:T1;
    in=(i-1)*7+1:(i-1)*7+7;
    ma(i,1)=std(msft(in,4));
    %the weekly return variance is equal to the sum of daily variances
    for j=1:7;
        t=(i-1)*7+j;
        ma(i,2)=ma(i,2)+(msft(t,4)^2)/7;
    end;
    ma(i,2) = sqrt(ma(i,2));
end;
disp(corrcoef(ma(:,1),ma(:,2)));

%plot the different time series
%weekly moving average
subplot(3,1,1); plot(ma(:,1));
%weekly sqrt of sum of daily squared returns
subplot(3,1,2); plot(ma(:,2));
%white noise
WN=exp(randn(T1,1)); subplot(3,1,3); plot(WN);

%% Predictability of Return Volatility

%predictability for white noise
X=[ones(T1-1,1) WN(1:T1-1,1)];
Y=WN(2:T1,1);
[regcoef sterr a3 a4 rsq]=regress(Y,X);
disp([regcoef sterr]);

%predictability for actual data
X=[ones(T1-1,1) ma(1:T1-1,2)];
Y=ma(2:T1,2);
[regcoef sterr a3 a4 rsq]=regress(Y,X);
disp([regcoef sterr]);
%Note regcoef(2) is positive and significant!
disp(rsq(1));
%Note, R2 is 30.39%!
%Even with a simple linear model we get
%predictability! Can you think of better models?

%scatter plots 
subplot(2,1,1); plot(WN(1:T1-1,1),WN(2:T1,1),'.'); 
subplot(2,1,2); plot(ma(1:T1-1,2),ma(2:T1,2),'.'); 

%% Simulation Stochastic Volatility

%%% Discuss Function simsecSV %%%

%simulating stochastic volatility
sigmaS=.0025; rho=.96;
r=simsecSV(mu,sigma,J1,J2,p1,p2,rho,sigmaS,T);
T1=floor(T/7); masim=zeros(T1,2);
for i=1:T1;
    in=(i-1)*7+1:(i-1)*7+7;
    masim(i,1)=std(r(in,1));
    for j=1:7;
       t=(i-1)*7+j;
       masim(i,2)=masim(i,2)+(r(t,1)^2)/7;
    end;
    masim(i,2)=sqrt(masim(i,2));
end;

%predictability regression and plotting the vol time series
X=[ones(T1-1,1) masim(1:T1-1,1)]; Y=masim(2:T1,1);
[regcoef sterr a3 a4 rsq]=regress(Y,X);
disp([regcoef sterr]); disp(rsq(1));
subplot(2,1,1); plot(ma(:,1)); title('msft 7 day vol');
subplot(2,1,2); plot(masim(:,1)); title('simulated 7 day vol');

