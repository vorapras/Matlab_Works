function C=BSCall(S,K,v,r,t)

d1=(log(S/K)+(r+v^2/2)*t)/v/t;
d2=d1-v*sqrt(t);
C=S*normcdf(d1)-exp(-r*t)*K*normcdf(d2);