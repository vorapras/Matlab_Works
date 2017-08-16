function out=pricederiv(pS,RfAgg,sigmaAgg,X,N)

sigma=sigmaAgg/sqrt(N);
Rf=RfAgg^(1/N);

%initialize price grids for the derivative (C) and underlying (P)
%N+1 is the time until expiry (time dimension)
%2^N is the number of possible different prices of underlying at expiry
C=zeros(2^N,N+1);
P=zeros(2^N,N+1);

%create a price grid for underlying
P(1,1)=pS;
for i=1:N;
    for j=1:2^(i-1);
        P((j-1)*2+1,i+1)=P(j,i)*(1+sigma);
        P((j-1)*2+2,i+1)=P(j,i)*(1-sigma);
        %disp([i j i+1 (j-1)*2+1 (j-1)*2+2]);
    end;
end;

%define payout at expiry of derivative as a function of underlying
C(:,N+1)=max(P(:,N+1)-X,0); %a european call
%C(:,N+1)=max(X-P(:,N+1),0); %a european put
%C(:,N+1)=abs(P(:,N+1)-X); %a volatility hedge
%we don't have to price only puts and calls, we can price
%any exotic option, as long as its payout depends on the 
%underlying stock

%create price grid for option
for k=1:N;
    i=N+1-k;
    for j=1:2^(i-1);
        Ch=C((j-1)*2+1,i+1);
        Cl=C((j-1)*2+2,i+1);
        pStemp=P(j,i);
        x=pricebinomial(pStemp,Rf,sigma,Ch,Cl);
        C(j,i)=x(1);
        %disp([i j i+1 (j-1)*2+1 (j-1)*2+2]);
    end;
end;
out=C(1,1);
        