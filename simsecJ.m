function r=simsecJ(mu,sigma,J1,J2,p1,p2,T);
x=randn(T,1); y=rand(T,1);
for t=1:T;
    if y(t)<p1; 
       J(t)=J1;
    elseif y(t)<p1+p2; 
       J(t)=J2;
    else
       J(t)=0;
    end;
    r(t,1)=exp(mu+sigma*x(t)+J(t))-1;
end;
