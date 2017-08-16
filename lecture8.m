%% Code for FM457A and FM457B - Lecture 8

% clear workspace
clear all;
close all;
clc;


%% Pricing a Call with No Arbitrage

%%% Discuss Function pricebinomial %%%
pS=100; Rf=1.02; sigma=.1; Ch=10; Cl=0;
pricebinomial(pS,Rf,sigma,Ch,Cl)

%% Binomial Trees

%%% Note: The assumption that the world only has two states is unrealistic.
%%% However, it's not untrealistic to assume that the price in one minute
%%% can only take on two values. In fact, at the limit, the binomial
%%% assumption implies a log-normal distribution of prices at expiry. Put
%%% differently, the binomial tree model converges to the Black&Scholes
%%% model in the limit (together with some other assumptions...).

%build a binomial tree in matrix form
N=3; P=zeros(2^N,N+1);
%create a price grid for underlying
P(1,1)=pS;
for i=1:N;
  for j=1:2^(i-1);
    P((j-1)*2+1,i+1)=P(j,i)*(1+sigma);
    P((j-1)*2+2,i+1)=P(j,i)*(1-sigma);
    disp([i j i+1 (j-1)*2+1 (j-1)*2+2]);
  end;
end;
disp(P);


%% Find the Call Price using Backward Induction

%Recursively define prices backwards
X=100; C(:,N+1)=max(P(:,N+1)-X,0); %call option 
for k=1:N;
    i=N+1-k;
    for j=1:2^(i-1);
        Ch=C((j-1)*2+1,i+1);
        Cl=C((j-1)*2+2,i+1);
        pStemp=P(j,i);
        out=pricebinomial(pStemp,Rf,sigma,Ch,Cl);
        C(j,i)=out(1);
        disp([i j i+1 (j-1)*2+1 (j-1)*2+2]);
    end;
end;
disp(C);

%%% Note: Call with longer time to maturity is more valuable

%% Visualize Convergence to BS Call Price

%%% Discuss Function pricederiv - i.e. the function consist out of the two
%%% steps introduced above: 1.Find tree with underlying's price 2.Find the
%%% European call price through backward induction

pS=100; Rf=1.02; sigmaAgg=.3; X=100;
for N=1:15;
    out(N,1)=N;
    out(N,2)=pricederiv(pS,Rf,sigmaAgg,X,N);
end;
plot(out(:,1),out(:,2));

blsprice(pS, X, Rf-1, 1, sigmaAgg)

%%% Note: the call price converges to BS call price as N grows

%% Investigating Call Prices for Different Strikes

%price and plot calls for different strikes
pS=100; Rf=1.02; sigmaAgg=.3; N=10;
for i=1:50;
    X=40+120*(i-1)/(50-1);   
    out(i,1)=X;
    out(i,2)=pricederiv(pS,Rf,sigmaAgg,X,N);
end;
plot(out(:,1),out(:,2));
xlabel('Strike'); ylabel('Call');

%%% Note: Call price is decreasing in the strike - makes sense since strike
%%% is what one needs to pay if option is exercised.

%% Investigating Call Prices for Different Prices of Underlying

%price and plot calls for different prices of underlying
X=100; Rf=1.02; sigmaAgg=.3; N=10;
for i=1:50;
    pS=40+120*(i-1)/(50-1);   
    out(i,1)=pS;
    out(i,2)=pricederiv(pS,Rf,sigmaAgg,X,N);
end;
plot(out(:,1),out(:,2));
xlabel('Price'); ylabel('Call');

%%% Note: Given a specific strike, the price of a call is increasing in the
%%% price of the underlying - makes sense, since this is what one gets if
%%% option is exercised.

%% Investigating Call Prices for Different Volatilities

%price and plot calls for different volatilities
X=100; Rf=1.02; pS=100; N=10;
for i=1:50;
    sigmaAgg=.01+.8*(i-1)/(50-1);   
    out(i,1)=sigmaAgg;
    out(i,2)=pricederiv(pS,Rf,sigmaAgg,X,N);
end;
plot(out(:,1),out(:,2));
xlabel('Sigma'); ylabel('Call');

%%% Note: Call price is increasing in volatility - makes sense since
%%% probability of the option ending up in the money increases with
%%% volatility.


