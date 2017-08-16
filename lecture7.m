%% Code for FM457A and FM457B - Lecture 7

% clear workspace
clear all;
close all;
clc;


%% Fama and French - 3 Factor Model

%load data
ff6ret; ff3fact;

%look at data
disp(mean(FF6ret(:,2:7))); % Returns are oredered as follows: Small companies with
% Low Middle High B/M columns 2:4; Large companies with Low Middle High B/M
% columns 5:7
disp(mean(FFfactors(:,2:5))); % Returns are ordered as follows Mkt-RF SMB HML RF

%% Test CAPM on FF6

%run regressions of the 6 portfolios on the market factor only and display
%the alphas - remember alpha = 0 needs to hold in order to justify CAPM.
[T w]=size(FFfactors);
X=[ones(T,1) FFfactors(:,2)];
for i=1:6;
    Y=FF6ret(:,1+i)-FFfactors(:,5);
    [regcoef sterr]=regress(Y,X);
    outcapm(i,:)=[regcoef(1) sterr(1,1) sterr(1,2)];
end;
disp(outcapm);
subplot;

%plot actual vs fitted returns. If CAPM is correct the returns should line
%up on the 45 degree line
plot(mean(FF6ret(:,2:7)),mean(FF6ret(:,2:7)),'k-'); hold on;
plot(mean(FF6ret(:,2:7)),mean(FF6ret(:,2:7))-outcapm(:,1)','b.');


%% %% Test 3 Factor Model on FF6

%run regressions of the 6 portfolios on the 3 factors and display
%the alphas - remember alpha = 0 needs to hold in order to justify CAPM.
[T w]=size(FFfactors);
X=[ones(T,1) FFfactors(:,2:4)];
for i=1:6;
    Y=FF6ret(:,1+i)-FFfactors(:,5);
    [regcoef sterr]=regress(Y,X);
    out3f(i,:)=[regcoef(1) sterr(1,1) sterr(1,2)];
end;
disp(out3f);

%plot actual vs fitted returns. If 3 model is correct the returns should line
%up on the 45 degree line
plot(mean(FF6ret(:,2:7)),mean(FF6ret(:,2:7)),'k-'); hold on;
plot(mean(FF6ret(:,2:7)),mean(FF6ret(:,2:7))-out3f(:,1)','rx');

%% Return Predictability

%load data - contains GE returns (incl. dividends) in the 3. column and
%returns (excl. dividends) in the 4. column
aggret;

%create a monthly time series for prices (Jan 1926 normalized to 1) and
%dividends
[T b]=size(aggreturn);
price(1)=1; dividend(1)=0;
for t=2:T;
    timeline(t)=round(aggreturn(t,2)/10000);
    price(t)=price(t-1)*(aggreturn(t,4)+1); 
    %calculate price of firm next year using the ex-dividend return
    dividend(t)=price(t-1)*(aggreturn(t,3)+1)-price(t); 
    %calculate dividend by subtracting ex-dividend price from 
    %cum-dividend price
end; 
dividend(1)=dividend(2); %use some number other than zero

%%
%load data - contains monthly data for Jan 1926 - Dec 2006. Its columns are
%date, price (normalized to Jan 1926 being 1), dividend, P/D, aggregate
%return
pdin;

% Gordon Growth: P/D high implies prices are high relative to payouts, so
% high growth is expected
T=length(pddata);
timeline=(pddata(1,1)+(pddata(T,1)-pddata(1,1))*([1:T]-1)/(T-1))/10000;
subplot(2,1,1); plot(pddata(:,1),pddata(:,5)); title('Aggregate Return'); 
in=1:12:T; 
subplot(2,1,2); plot(pddata(in,1),pddata(in,4),'LineWidth',2);
title('Price/Dividend');

%%% Note: P/D has increased over the period. Have distributions changed so
%%% dramatically? 
%%% P/E is a better measure, since investors are getting paid in ways other
%%% than dividends.

%% Monthly Predicability Regressions

[T a]=size(pddata);
Y=pddata(7:T,5); 
X=[ones(T-6,1) pddata(1:T-6,4)/12];
[regcoef sterr a3 a4 rsq]=regress(Y,X);
disp([regcoef sterr]); disp(rsq(1));

%%% Note: Regression coefficient is negative (high preices today lead to
%%% low returns in the future), but not significant, R^2 is very low

%% Annual Predicability Regressions
T1=floor(T/12); %number of months in our data
for i=2:T1;
    in=(i-1)*12+1:(i-1)*12+12; %selective index for year i
    Y(i,1)=sum(pddata(in,5));
    X(i,1)=1;
    X(i,2)=pddata((i-1)*12+1-7,4)/12;
end;
[regcoef sterr a3 a4 rsq]=regress(Y(2:T1),X(2:T1,:));
disp([regcoef sterr]); disp(rsq(1));

%%% Note: Regression coefficient is negative (high preices today lead to
%%% low returns in the future) and R^2 is 7.3%. Annual returns over this
%%% sample period are predicted by P/D!

