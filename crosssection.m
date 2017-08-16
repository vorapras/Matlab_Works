function R=crosssection(mu,sigma,rf,B,sigmai,T)
[a nsec]=size(B);
R=zeros(T,nsec); 
X=randn(T,1); Xi=randn(T,nsec);
rm=mu+sigma*X;
mui=rf+B*(mu-rf);
for i=1:nsec;
    R(:,i)=mui(i)+(rm-mu).*B(i)+Xi(:,i)*sigmai(i);
end;

% mu: expected return on the market
% sigma: colatility of the market
% rf: risk-free rate
% mui: expected return of the risky assets
% B: vector containing the betas of the risky assets
% sigmai: vector containing the volatilities of the risky assets
% T: number of time periods
