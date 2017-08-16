function out=pricebinomial(pS,Rf,sigma,Ch,Cl);

D=(Ch-Cl)/(pS*2*sigma);
B=(Ch-(1+sigma)*pS*D)/Rf;
pC=B+pS*D;
out=[pC D];