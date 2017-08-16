function r=simsecSV(mu,muS,J1,J2,p1,p2,rho,sigmaS,T);
x=randn(T,1); y=rand(T,1); z=randn(T+1,1); W(1)=0; 
J=zeros(T,1); r=zeros(T,1);
for t=1:T;
    W(t+1)=rho*W(t)+sigmaS*z(t+1);
    sigma(t)=abs(muS+W(t));
    if y(t)<p1; 
       J(t)=J1;
    elseif y(t)<p1+p2; 
       J(t)=J2;
    else
       J(t)=0;
    end;
    r(t,1)=exp(mu+sigma(t)*x(t)+J(t))-1;
end;