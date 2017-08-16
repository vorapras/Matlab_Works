function r=simsec(mu,sigma,T);
x=randn(T,1); 
for t=1:T;
    r(t,1)=exp(mu+sigma*x(t))-1;
end;
