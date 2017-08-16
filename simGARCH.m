%this function simulates GARCH(1,1) w/ params w,a,b
%output is Tx2 matrix with r in 1st column, sigmasq in 2nd
function out=simGARCH(w,a,b,T);
r=zeros(T+100,1); 
sigmasq=zeros(T+100,1);  
r(1)=0; sigmasq(1)=w/(1-a-b); 
for t=1:T-1+100; 
    sigmasq(t+1)=w+a*(r(t)^2)+b*sigmasq(t); 
    r(t+1)=randn(1)*sqrt(sigmasq(t+1)); 
end;
out=[r(101:T+100) sigmasq(101:T+100)];

%this function incorporates a burning period of 100 periods.
