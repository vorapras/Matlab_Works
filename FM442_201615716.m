% Multivariate GARCH. Select two return series. Estimate a univariate GARCH model of each
% return, and a multivariate GARCH model of both returns, and analyze the results. Focus on
% the statistical significance between the univariate and multivariate models, the difference in
% the univariate and multivariate volatility forecasts and residuals, and the behavior of the
% covariance over time.
%% *** DataProcess ***
clc;
clear all;
clf;
close all;

%% *** Part 1: Loading daily data ***
Data_D           = importdata('FM442_201615716.csv');
Tickers          = ['MSFT'; 'XOM '; 'GE  '; 'JPM '; 'INTC'; 'C   '; ...
                       'SPX '];
%% *** Part 2 - Cleaning up and organizing the data ***
% Retrieving information about data
[NRows, NCols]   = size(Data_D.data);

% Preallocating space for matrices with PermNo and Dates
PermNo           = NaN(NRows, 1);
Dates            = NaN(NRows, 1);

for i = 1:NRows
    % Note: the str2num function converts text to data
    PermNo(i,1)  = str2num(Data_D.textdata{i+1, 1});
    Dates( i,1)  = str2num(Data_D.textdata{i+1, 2});
end

% Collecting data in a single matrix
DataNum          = [Dates, PermNo, Data_D.data];

% Getting unique values of stock identifiers and dates
PermNoUnique     = unique(PermNo);
DatesUnique      = unique(Dates);

% Computing number of stocks and dates
NStocks          = length(PermNoUnique);
NDates           = length(DatesUnique);

%% *** Part 3 - Computing total returns and total returns index ***

% Allocate memory for total returns for each stock
TR_Stocks        = NaN(NDates, NStocks);

for i = 1:NStocks
   % Extract the rows that correspond to that stock
   Indices       = find(PermNo==PermNoUnique(i));
   StockData     = DataNum(Indices, :);
   
   % Sort the data according to date
   StockData     = sortrows(StockData);

   % Store total returns for this stock in matrix
   TR_Stocks(:, i) = StockData(:, 3);
   
   % Store total returns for S&P 500
   if i == 1
       TR_SPX      = StockData(:, 4);
   end
end

% Compiling single matrix with all total returns data
TR                 = [TR_Stocks, TR_SPX];  

% Computing total returns indices
TRIndices          = NaN(NDates + 1, NStocks + 1);
TRIndices(1, :)    = ones(1, NStocks + 1);

for i=1:NDates
    TRIndices(i+1, :) = TRIndices(i, :) .* (1 + TR(i, :));
end

% Computing log returns
LogReturns         = log(1 + TR);

%% *** Part 4 - Divide data into 2 period  ***
% Estimation Period 2005 to 2007 and Forecasting Period 2008

% Burn-in Period
Window              = 253;
Indices             = [1 2 6]; % MSFT XOM C

% Set up Estimation data from 20050103 to 20081231
MyData              = LogReturns(3785:4791,Indices);
MyDataSqr           = MyData.*MyData ;
PortData            = MyData*[1/3;1/3;1/3];
EstData             = [MyData,PortData];

% Find Length of All Data
DataLength          = length(EstData);

%% *** Estimating GARCH and its family model ***

%% *** Part 5 - Estimating normal GARCH from 20050103 to 20081231 ***
% Initialise Data
GARCHData           = EstData  ;
GARCHDataSq         = GARCHData .* GARCHData;
% Initialise space for estimation
GARCHCoeffs         = NaN(4, 4);
GARCHModel          = garch(1, 1);
GARCHlogL           = NaN(1,4);
GARCHAIC            = NaN(1,4);
GARCHBIC            = NaN(1,4);

for i=1:4
% Estimate GARCH model and Log-Likelihood
[GARCHPara,~,GARCHlogL(1,i),~] = estimate(GARCHModel, GARCHData(:,i));
% Calculate AIC and BIC
[aic,bic]           = aicbic(GARCHlogL(1,i), 2, DataLength);
GARCHAIC(1,i)       = aic ;
GARCHBIC(1,i)       = bic ;
end

% Extract coefficients
GARCHCoeffs(1,:)    = GARCHPara.Constant;
GARCHCoeffs(2,:)    = GARCHPara.ARCH{1};
GARCHCoeffs(3,:)    = GARCHPara.GARCH{1};
GARCHCoeffs(4,:)    = GARCHPara.UnconditionalVariance;

% Compute conditional volatility estimates
GARCHVar            = NaN(DataLength, 4);
GARCHVar(1, :)      = GARCHCoeffs(4,:) ;

for t=2:DataLength
   GARCHVar(t, :)   = GARCHCoeffs(1, :) ...
                         + GARCHCoeffs(2,i) .* GARCHDataSq(t-1, :) ...
                         + GARCHCoeffs(3,i) .* GARCHVar(t-1, :);
end

% Not include burned-in period
GARCHVol            = sqrt(GARCHVar((Window+1):end, :));
% Extract each stock conditional volatility
MSFTGARCHVol        = GARCHVol(:, 1);
XOMGARCHVol         = GARCHVol(:, 2);
CGARCHVol           = GARCHVol(:, 3);
PortRetGARCHVol     = GARCHVol(:, 4);

% Calculate Residual from Portfolio Return
GARCHResid          = PortData(Window+1,1)./PortRetGARCHVol;

% Diagnose distribution of Residual by QQ plots
% QQ plots versus Normal Distribution
a=1; % first graph
figure(a);
qqplot(GARCHResid)
a= a+1;

% Calculate MSE to compare forecasting accuracy
PortRetSq           = PortData.*PortData;

for i = 1:DataLength   
    Diff            = 0;
    Diff            = Diff+(GARCHVar(i,4)- PortRetSq(i,1))^2 ;
end
GARCHMSE            = Diff/DataLength ;

%% *** Part 6  - Estimating GJR-GARCH from 20050103 to 20081231***
% Initialise space for estimation
GJRGARCHCoeffs      = NaN(4, 4);
GJRGARCHModel       = gjr(1, 1);
GJRGARCHlogL        = NaN(1,4);
GJRGARCHAIC         = NaN(1,4);
GJRGARCHBIC         = NaN(1,4);


for i=1:4
% Note: parameters are, in order, omega, alpha, gamma and beta
[GJRGARCHPara,~,GJRGARCHlogL(1,i),~] = estimate(GJRGARCHModel,GARCHData(:, i));
% Extracting parameters
GJRGARCHCoeffs(1,i)    = GJRGARCHPara.Constant;
GJRGARCHCoeffs(2,i)    = GJRGARCHPara.ARCH{1};
GJRGARCHCoeffs(3,i)    = GJRGARCHPara.GARCH{1};
GJRGARCHCoeffs(4,i)    = GJRGARCHPara.Leverage{1};
% Calculate AIC and BIC
[aic,bic]        = aicbic(GJRGARCHlogL(1,i), 3, DataLength);
GJRGARCHAIC(1,i) = aic ;
GJRGARCHBIC(1,i) = bic ;
end

% Computing volatility estimates for each series
% First, allocate memory for variance estimates
GJRGARCHVar      = NaN(DataLength, 4);

% Will use the first year as a burn-in period, so volatility estimates
% are not sensitive to how we estimate the first-period variance

% Can use squared returns or GARCH-like estimates (our choice)
GJRGARCHVar(1, :)= GJRGARCHCoeffs(1,:) ./ (ones(1, 4) - GJRGARCHCoeffs(2,:) - GJRGARCHCoeffs(4,:));

% Creating matrix with the sign of returns
ReturnsIsNeg     = ones(size(GARCHData));
ReturnsIsNeg(GARCHData > 0) = 0;

% Loop to compute variance estimates
for t=2:DataLength
    GJRGARCHVar(t, :) = GJRGARCHVar(1, :) + (GJRGARCHCoeffs(2,:) + GJRGARCHCoeffs(3,:) .* ReturnsIsNeg(t-1, :)) ...
                        .* GARCHDataSq(t-1, :) + GJRGARCHCoeffs(4,:).* ...
                        GJRGARCHVar(t-1, :);
end

% Compute volatility, discard first 253 observations
GJRGARCHVol            = sqrt(GJRGARCHVar);
GJRGARCHVol(1:253, :)  = [];

% Extract each stock conditional volatility
MSFTGJRGARCHVol        = GJRGARCHVol(:, 1);
XOMGJRGARCHVol         = GJRGARCHVol(:, 2);
CGJRGARCHVol           = GJRGARCHVol(:, 3);
PortRetGJRGARCHVol     = GJRGARCHVol(:, 4);

% Calculate Residual from Portfolio Return
GJRGARCHResid          = PortData(Window+1,1)./PortRetGJRGARCHVol;

% Diagnose distribution of Residual by QQ plots
% QQ plots versus Normal Distribution
figure(a);
qqplot(GJRGARCHResid)
a = a+1;

% Calculate MSE to compare forecasting accuracy
for i = 1:DataLength   
    Diff            = 0;
    Diff            = Diff+(GJRGARCHVar(i,4)- PortRetSq(i,1))^2 ;
end
GJRGARCHMSE         = Diff/DataLength ;

%% *** Part 7  - Estimating normal EGARCH from 20050103 to 20081231 ***
EGARCHCoeffs           = NaN(4, 4);
EGARCHModel            = egarch(1, 1);
EGARCHlogL             = NaN(1,4);
EGARCHAIC              = NaN(1,4);
EGARCHBIC              = NaN(1,4);

for i=1:4
 % Note: parameters are, in order, omega, alpha, gamma and beta
 [EGARCHPara,~,EGARCHlogL(1,i),~] = estimate(EGARCHModel,GARCHData(:, i));
 % Extracting parameters
 EGARCHCoeffs(1,i)    = EGARCHPara.Constant;
 EGARCHCoeffs(2,i)    = EGARCHPara.ARCH{1};
 EGARCHCoeffs(3,i)    = EGARCHPara.GARCH{1};
 EGARCHCoeffs(4,i)    = EGARCHPara.Leverage{1};
 % Calculate AIC and BIC
 [aic,bic]           = aicbic(EGARCHlogL(1,i), 3, DataLength);
 EGARCHAIC(1,i)      = aic ;
 EGARCHBIC(1,i)      = bic ;
end

% Computing volatility estimates for each series
% First, allocate memory for variance estimates
logEGARCHVar         = NaN(DataLength, 4);
Epsilon              = NaN(DataLength-1, 4);
ExpectedEp           = sqrt(2/pi);
% Can use squared returns or GARCH-like estimates (our choice)
logEGARCHVar(1, :)   = EGARCHCoeffs(1,:) ./ (ones(1, 4) - GJRGARCHCoeffs(2,:) - GJRGARCHCoeffs(3,:));

% Loop to compute variance estimates
for t=2:DataLength
    Epsilon(t-1,:)    = GARCHData(t-1,:)/sqrt(exp(logEGARCHVar(t-1, :)));
    logEGARCHVar(t, :) = EGARCHCoeffs(1,:) + EGARCHCoeffs(4,:).* Epsilon(t-1,:)...
                        + EGARCHCoeffs(3,:).*logEGARCHVar(t-1,:)...
                        + EGARCHCoeffs(2,:).*(abs(Epsilon(t-1,:))-ExpectedEp); 
end

% Compute volatility, discard first 253 observations
EGARCHVar           = exp(logEGARCHVar);
EGARCHVol           = sqrt(EGARCHVar);
EGARCHVol(1:253, :) = [];

% Extract each stock conditional volatility
MSFTEGARCHVol       = EGARCHVol(:, 1);
XOMEGARCHVol        = EGARCHVol(:, 2);
CEGARCHVol          = EGARCHVol(:, 3);
PortRetEGARCHVol    = EGARCHVol(:, 4);

% Calculate Residual from Portfolio Return
EGARCHResid          = PortData(Window+1,1)./PortRetEGARCHVol;

% Diagnose distribution of Residual by QQ plots
% QQ plots versus Normal Distribution
figure(a);
qqplot(EGARCHResid)
a = a+1;

% Calculate MSE to compare forecasting accuracy
for i = 1:DataLength   
    Diff            = 0;
    Diff            = Diff+(EGARCHVar(i,4)- PortRetSq(i,1))^2 ;
end
EGARCHMSE         = Diff/DataLength ;

%% *** Part 8  - Estimating Correlation from the actual data ***
% Computing rolling 20-day correlations
MultCorr            = NaN(DataLength, 3);
for t=(Window+1): DataLength  
    % Gathering data rolling 20-day for each stock return
    Data1           = MyData((t-20):(t-1), 1); % MSFT 
    Data2           = MyData((t-20):(t-1), 2); % XOM 
    Data3           = MyData((t-20):(t-1), 3); % C
    
    CorrMatMSFT_XOM = corrcoef(Data1, Data2);
    CorrMatMSFT_C   = corrcoef(Data1, Data3);
    CorrMatXOM_C    = corrcoef(Data2, Data3);
   
    MultCorr(t, 1)  = CorrMatMSFT_XOM(1,2);
    MultCorr(t, 2)  = CorrMatMSFT_C(1,2);
    MultCorr(t, 3)  = CorrMatXOM_C(1,2);    
end
% Not include Burn-in-period
MultCorr(1:Window,:)= [];

%% *** Part 9 - Calculate Value-at-Risk  ***
PortCondVol         = PortRetGARCHVol ; %PortRetGARCHVol ; % PortRetGJRGARCHVol ; % PortRetEGARCHVol ;
p                   = 0.01;            % VaR level
InvCDF              = norminv(p);
ValueAtRiskEstG     = - InvCDF * PortCondVol;

% Set up Violation Ratios function
Breach              = (PortData((Window+1):end) < - ValueAtRiskEstG);
LenBreach           = size(Breach, 1);
pHat                = sum(Breach)/size(Breach, 1);
VRtest              = pHat/p;

% Unconditional coverage test
V1                  = sum(Breach);
V0                  = size(Breach, 1) - V1;
LogLConst           = V1 * log(p) + V0 * log(1-p);
LogLUnconst         = V1 * log(pHat) + V0 * log(1-pHat);
TestStatisticUC     = -2 * (LogLConst - LogLUnconst);
SignificanceUC      = 1 - chi2cdf(TestStatisticUC, 1);

% Conditional coverage test
% Estimating the transition probabilities
Aux01               = (ones(LenBreach - 1, 1) - Breach(1:(end -1), :)) ...
                         .* Breach(2:end, :);
p01                 = sum(Aux01) / (LenBreach - 1 ...
                         - sum(Breach(1:(end -1), :)));
                     
Aux11               = Breach(1:(end -1), :) .* Breach(2:end, :);
p11                 = sum(Aux11) / sum(Breach(1:(end -1), :));

% Computing log-likelihood; note that the constrained version here is 
% equal to the unconstrained version above
LogLConst           = LogLUnconst;

% Running loop to compute the unconstrained log-likelihood
LogLUnconst         = Breach(1, 1) * log(pHat) + ...
                         (1 - Breach(1, 1)) * log(1 - pHat);
for i=2:size(Breach, 1)
    if     (Breach(i-1, 1) == 0 && Breach(i, 1) == 0)
        LogLUnconst = LogLUnconst + log(1 - p01);
    elseif (Breach(i-1, 1) == 0 && Breach(i, 1) == 1)
        LogLUnconst = LogLUnconst + log(p01);
    elseif (Breach(i-1, 1) == 1 && Breach(i, 1) == 0)
        LogLUnconst = LogLUnconst + log(1 - p11);
    elseif (Breach(i-1, 1) == 1 && Breach(i, 1) == 1)
        LogLUnconst = LogLUnconst + log(p11);
    end
end

% Chi-Square Test
TestStatisticCC     = -2 * (LogLConst - LogLUnconst);
SignificanceCC      = 1 - chi2cdf(TestStatisticCC, 1);

% Calculate VaR at 1%,5% and 10%
alpha               = [0.01,0.05,0.1]; % VaR level
InvCDF              = norminv(alpha);
VaRGARCH            = NaN(size(PortCondVol,1),3);
for i = 1:3
VaRGARCH(:,i)       = - InvCDF(1,i) * PortCondVol;
end

%% *** Estimating Multivariate GARCH and its family model from 20050103 to 20081231 ***

%% *** Part 10: BEKK with stock data ***
% Estimate BEKK(1, 1) - Do not run this since it takes so long
[BEKKParms,BEKKlogL,BEKKVarMat,~]      = bekk(MyData, [], 1, 0, 1, 'Full');

% Section to save and retrieve BEKKParms, if needed
save('BEKKParms.mat', 'BEKKParms');
save('BEKKlogL.mat', 'BEKKlogL');
load('BEKKParms');
load('BEKKlogL');

% Calculate AIC and BIC
BEKKaic            = (2*24)-(2*BEKKlogL);
BEKKbic            = (log(DataLength)*24)-(2*BEKKlogL);

% Convert parameters to usual format
[BEKK_C, BEKK_A, ~, BEKK_B] = bekk_parameter_transform(BEKKParms, ...
                      1, 0, 1, 3, 3);

% Compute variance matrix estimates at each point in time -
% use first 253 observations to initialize computation.

% The following matrix will hold the estimates; first row is for the
% initial period
BEKKVarMat         = NaN(3, 3, DataLength);
BEKKVarMat(:,:, 1) = cov(MyData(1:Window, :));

for i=2:DataLength
    LogRet            = MyData(i-1, :);
    BEKKVarMat(:,:,i) = BEKK_C ...
          + BEKK_A' * (LogRet' * LogRet) * BEKK_A ...
          + BEKK_B' * BEKKVarMat(:,:,i-1) * BEKK_B;
end

% Delete first year of estimates
BEKKVarMat         = BEKKVarMat(:, :, (Window+1):end);

% Extract variance estimates for each stock, as well as correlations
BEKKVolMSFT        = sqrt(BEKKVarMat(1, 1, :));
BEKKVolXOM         = sqrt(BEKKVarMat(2, 2, :));
BEKKVolC           = sqrt(BEKKVarMat(3, 3, :));

% Extracting correlation information, deleting first year
% MSFT XOM 
BEKKCorr12           = BEKKVarMat(1, 2, :) ./ (BEKKVolMSFT .* BEKKVolXOM);
% MSFT C
BEKKCorr13           = BEKKVarMat(1, 3, :) ./ (BEKKVolMSFT .* BEKKVolC);
% XOM C
BEKKCorr23           = BEKKVarMat(2, 3, :) ./ (BEKKVolXOM .* BEKKVolC);

% Reshape matrices into vectors
BEKKVolMSFT         = reshape(BEKKVolMSFT, DataLength - Window, 1);
BEKKVolXOM          = reshape(BEKKVolXOM, DataLength - Window, 1);
BEKKVolC            = reshape(BEKKVolC, DataLength - Window, 1);
BEKKCorr12          = reshape(BEKKCorr12, DataLength - Window, 1);
BEKKCorr13          = reshape(BEKKCorr13, DataLength - Window, 1);
BEKKCorr23          = reshape(BEKKCorr23, DataLength - Window, 1);

% Produce charts comparing output of BEKK with GARCH information
% Creating time labels
StartDate             = 2006;
EndDate               = 2009; 
TimeLabels            = linspace(StartDate, EndDate, DataLength - Window);

% Produce figure
figure(a);
subplot(3, 1, 1);
titleStr = 'Volatility Estimates for MSFT';
title(titleStr);
xlabel('Dates');
axis([StartDate EndDate 0 1]);
grid on;
hold on;
plot(TimeLabels, sqrt(252) * MSFTGARCHVol, 'LineStyle', '-' , ...
                  'LineWidth', 1, 'Color', 'blue');
plot(TimeLabels, sqrt(252) * BEKKVolMSFT, 'LineStyle', '--', ...
                 'LineWidth', 1, 'Color', 'red');
legend('GARCH', 'BEKK', 'location','best');              
set(gcf, 'color', 'white'); 
ylabel('Volatility');

subplot(3, 1, 2);
titleStr = 'Volatility Estimates for EXXON';
title(titleStr);
xlabel('Dates');
axis([StartDate EndDate 0 1.5]);
grid on;
hold on;
plot(TimeLabels, sqrt(252) * XOMGARCHVol, 'LineStyle', '-' , ...
                  'LineWidth', 1, 'Color', 'blue');
plot(TimeLabels, sqrt(252) * BEKKVolXOM , 'LineStyle', '--', ...
                 'LineWidth', 1, 'Color', 'red');
legend('GARCH', 'BEKK', 'location','best');              
set(gcf, 'color', 'white'); 
ylabel('Volatility');

subplot(3, 1, 3);
titleStr = 'Volatility Estimates for C';
title(titleStr);
xlabel('Dates');
axis([StartDate EndDate 0 1.5]);
grid on;
hold on;
plot(TimeLabels, sqrt(252) * CGARCHVol, 'LineStyle', '-' , ...
                  'LineWidth', 1, 'Color', 'blue');
plot(TimeLabels, sqrt(252) * BEKKVolC, 'LineStyle', '--', ...
                 'LineWidth', 1, 'Color', 'red');
legend('GARCH', 'BEKK', 'location','best');              
set(gcf, 'color', 'white'); 
ylabel('Volatility');
a=a+1;

% Correlation Plot
% Produce figure
figure(a);
% MSFT XOM 
subplot(3, 1, 1);
titleStr = 'Correlation Estimates for  MSFT and XOM  ';
title(titleStr);
xlabel('Dates');
axis([StartDate EndDate -1.2 1.2]);
grid on;
hold on;
plot(TimeLabels, MultCorr(:,1), 'LineStyle', '-' , ...
                  'LineWidth', 1, 'Color', 'blue');
plot(TimeLabels, BEKKCorr12, 'LineStyle', '--', ...
                 'LineWidth', 1, 'Color', 'red');
legend('20-Day', 'Out-of-Sample DCC', 'location','best');              
set(gcf, 'color', 'white'); 
ylabel('Correlation');

% MSFT C
subplot(3, 1, 2);
titleStr = 'Correlation Estimates for  MSFT and C  ';
title(titleStr);
xlabel('Dates');
axis([StartDate EndDate -1.2 1.2]);
grid on;
hold on;
plot(TimeLabels, MultCorr(:,2), 'LineStyle', '-' , ...
                  'LineWidth', 1, 'Color', 'blue');
plot(TimeLabels, BEKKCorr13, 'LineStyle', '--', ...
                 'LineWidth', 1, 'Color', 'red');
legend('20-Day', 'Out-of-Sample DCC', 'location','best');              
set(gcf, 'color', 'white'); 
ylabel('Correlation');

% XOM C
subplot(3, 1, 3);
titleStr = 'Correlation Estimates for  XOM and C  ';
title(titleStr);
xlabel('Dates');
axis([StartDate EndDate -1.2 1.2]);
grid on;
hold on;
plot(TimeLabels, MultCorr(:,3), 'LineStyle', '-' , ...
                  'LineWidth', 1, 'Color', 'blue');
plot(TimeLabels, BEKKCorr23, 'LineStyle', '--', ...
                 'LineWidth', 1, 'Color', 'red');
legend('20-Day', 'Out-of-Sample DCC', 'location','best');              
set(gcf, 'color', 'white'); 
ylabel('Correlation');

a=a+1;

%% *** Part 11 - Computing portfolio volatility ***
Weights             = [1/3 ; 1/3 ; 1/3];
PortLogRet          = PortData ;
PortRetSq           = PortData.*PortData ;

% Computing conditional portfolio variance
PortCondVarBEKK      = NaN(DataLength - Window, 1);
for i=1: (DataLength - Window)
   CurrentVarMat    = BEKKVarMat(:, :, i);
   PortCondVarBEKK(i,:) = Weights' * CurrentVarMat * Weights;
end

% Computing conditional standard deviation, dropping initial window
ConditionalVol      = sqrt(PortCondVarBEKK);

% Calculate MSE to compare forecasting accuracy 
for i = 1:DataLength-Window
    Diff             = 0;
    Diff             = Diff+(PortCondVarBEKK(i,1)- PortRetSq(i+253,1))^2 ;
end
BEKKMSE               = Diff/DataLength ;

%% *** Part 12: CCC with stock data ***
% Estimate CCC model
CCCEstData        = MyData ;
CCCEstDataSq      = MyData.* MyData;

[CCCParms,CCClogL,CCCVarMat,~] = ccc_mvgarch(MyData, [], 1, 0, 1);

CCCOmega         = [CCCParms(1) CCCParms(4) CCCParms(7)];
CCCAlpha         = [CCCParms(2) CCCParms(5) CCCParms(8)];
CCCBeta          = [CCCParms(3) CCCParms(6) CCCParms(9)];
CCCRho           = [CCCParms(10) CCCParms(11) CCCParms(12)];

% Calculate AIC and BIC
CCCaic            = (2*12)-(2*CCClogL);
CCCbic            = (log(DataLength)*12)-(2*CCClogL);

% Constructing volatility estimates from GARCH Parameters
CCCVar           = NaN(DataLength, 3);
CCCVar(1, :)     = CCCOmega ./ (1 - CCCAlpha - CCCBeta);

for i=2:DataLength
    CCCVar(i, :) = CCCOmega + ...
                         CCCAlpha .* CCCEstDataSq(i - 1, :) + ...
                         CCCBeta .* GARCHVar(i - 1, 1 : 3);
end

% Convert to volatility, drop first year
CCCVol           = sqrt(CCCVar((Window+1):end,:));
% Putting together information for chart
CCCVolMSFT         = CCCVol(:,1);
CCCVolXOM          = CCCVol(:,2);
CCCVolC            = CCCVol(:,3);

% Calculate Epsilon from actual data
CCCEstEps          = CCCEstData((Window+1):end,:)./ CCCVol;

% Diagnose distribution of Residual by QQ plots
% QQ plots versus Normal Distribution
figure(a);
qqplot(CCCEstEps(:,1))
a = a+1;
figure(a);
qqplot(CCCEstEps(:,2))
a = a+1;
figure(a);
qqplot(CCCEstEps(:,3))
a = a+1;

%% *** Part 13 - Computing portfolio volatility ***
% Computing conditional portfolio variance
PortCondVarCCC       = NaN(DataLength - Window, 1);
for i=1: (DataLength - Window)
   CurrentVarMat    = CCCVarMat(:, :, i);
   PortCondVarCCC (i,:) = Weights' * CurrentVarMat * Weights;
end

% Computing conditional standard deviation, dropping initial window
ConditionalVol      = sqrt(PortCondVarCCC);

% Calculate MSE to compare forecasting accuracy 
for i = 1:DataLength-Window
    Diff             = 0;
    Diff             = Diff+(PortCondVarCCC(i,1)- PortRetSq(i+253,1))^2 ;
end
CCCMSE               = Diff/DataLength ;


%% *** Part 14 - Estimating parameters DCC  ***
% Use data from 20050103 to 20081231 
DCCEstData        = MyData ;
DCCEstDataSq      = DCCEstData .* DCCEstData;
NStocks           =  3;

% Estimate DCC model
[DCCParms,DCClogL,~,~] = dcc(DCCEstData, [], 1, 0, 1);

% Calculate AIC and BIC
DCCaic            = (2*14)-(2*DCClogL);
DCCbic            = (log(DataLength)*14)-(2*DCClogL);

% Preallocating space for storing DCC Parameters
DCCPar           = NaN(3,NStocks);

% Define the indices to read the value from DCC parameters 
for i=1:3
    DCCIndice    = i:3:i+3*(NStocks-1);
    DCCPar(i,:)  = DCCParms(DCCIndice);
end

% Store all values of each stocks to calculate variance   
DCCOmega         = DCCPar(1,:);
DCCAlpha         = DCCPar(2,:);
DCCBeta          = DCCPar(3,:);
DCCRho           = DCCParms(:,3*NStocks+1:3*NStocks+NStocks*(NStocks-1)/2);
DCC_a            = DCCParms(end-1);
DCC_b            = DCCParms(end);

% Converting data in DCCRho(Correlation) from vector to matrix 3x3
% Define the length of Rho Matrix
dim              = length(DCCRho);
RhoDim           = (-1+sqrt(1+8*dim))/2 ;
DCCRhoMat        = zeros(RhoDim ,RhoDim);

% Fill data from DCCRho in the right position
for k=1:RhoDim 
    for j=k:RhoDim 
        index           = sum(RhoDim :-1:RhoDim -k+2)+j-k+1;  
        DCCRhoMat(k,j)  = DCCRho(index);       
    end
end
 
% Transform matrix to find DCCRBar
DCCRhoMat2       = DCCRhoMat'; % transpose
% Add zeros row below
DCCRhoMat        = [DCCRhoMat; zeros(1,size(DCCRhoMat,1))];
% Add zeros column ahead
DCCRhoMat        = [zeros(size(DCCRhoMat,1),1) DCCRhoMat];
% Add zeros row up
DCCRhoMat2       = [zeros(1,size(DCCRhoMat2,1));DCCRhoMat2];
% Add zeros last column
DCCRhoMat2       = [DCCRhoMat2 zeros(size(DCCRhoMat2,1),1)];
% Create diagonal matrix
DiagMat          = diag(ones(NStocks,1));
% Summation all the matrix
DCCRBar          = DCCRhoMat+DCCRhoMat2+DiagMat;

%% *** Part 15 - Computing DCC Var-Cov Matrix in estimation period from 20050103 to 20081231***

% Preallocating DCC var-covariance matrix
DCCVarEstMAt     = NaN(3, 3, DataLength-Window); % Estimation Period

% Computing GARCH volatilities from DCC Parameters in 20050103 to 20071231 (or 20090102 to 20141231) 
% for each asset and standardized residuals
DCCEstVar        = NaN(DataLength, 3);
DCCEstVar(1, :)  = DCCOmega ./ (1 - DCCAlpha - DCCBeta);

for t=2:DataLength
 DCCEstVar(t, :) = DCCOmega + ...
                  DCCAlpha .* DCCEstDataSq(t - 1, :) + ...
                  DCCBeta .* DCCEstVar(t - 1, :);
end

% Not include burned-in period
DCCEstVol        = sqrt(DCCEstVar((Window+1):end,:));
% Calculate Epsilon from actual data
DCCEstEps        = DCCEstData((Window+1):end,:)./ DCCEstVol;

% Diagnose distribution of Residual by QQ plots
% QQ plots versus Normal Distribution
figure(a);
qqplot(CCCEstEps(:,1))
a = a+1;
figure(a);
qqplot(CCCEstEps(:,2))
a = a+1;
figure(a);
qqplot(CCCEstEps(:,3))
a = a+1;

% Compute value for DCC Qt and Rt ; 
% Initialize with estimate for RBar
DCCQt         = NaN(3, 3, DataLength-Window);
DCCRt         = NaN(3, 3, DataLength-Window);
DCCQt(:, :, 1)= DCCRBar;
DCCRt(:, :, 1)= DCCRBar;

% Calculate DCC covaraince matrix for each period to be further used in VaR
for t=2:DataLength-Window
   DCCQt(:, :, t) = (1 - DCC_a - DCC_b) * DCCRBar ...
                    + DCC_a * DCCEstEps(t-1, :)' * DCCEstEps(t-1, :) ...
                    + DCC_b * DCCQt(:, :, t-1);
   AuxEstMat         = [sqrt(DCCQt(1, 1, t))  sqrt(DCCQt(2, 2, t)) sqrt(DCCQt(3, 3, t))];
   % Correlation Matrix
   DCCRt(:, :, t) = DCCQt(:, :, t) ./ (AuxEstMat' * AuxEstMat);
   % Diagonal Varince Matrix
   DEstMat           = [DCCEstVol(t, 1) 0 0; 0 DCCEstVol(t, 2) 0; 0 0 DCCEstVol(t, 3)];
   % Covariance Matrix
   DCCVarEstMAt(:, :, t) = DEstMat * DCCRt(:, :, t) * DEstMat;
end

%% *** Part 16 - Compute historical volatilities and correlations ***

% Extracting correlation information, deleting first year
% MSFT XOM 
DCCCorr12          = DCCRt(1, 2, :);
DCCCorr12          = reshape(DCCCorr12,DataLength-Window,1);
% MSFT C
DCCCorr13          = DCCRt(1, 3, :);
DCCCorr13          = reshape(DCCCorr13, DataLength-Window, 1);
% XOM C
DCCCorr23          = DCCRt(2, 3, :);
DCCCorr23          = reshape(DCCCorr23, DataLength-Window, 1);

% Putting together information for chart
DCCVolMSFT         = DCCEstVol(:,1);
DCCVolXOM          = DCCEstVol(:,2);
DCCVolC            = DCCEstVol(:,3);

%% *** Part 17 - Computing portfolio volatility ***
% Computing conditional portfolio variance
PortCondVarDCC       = NaN(DataLength - Window, 1);
for i=1: (DataLength - Window)
   CurrentVarMat    = DCCVarEstMAt(:, :, i);
   PortCondVarDCC (i,:) = Weights' * CurrentVarMat * Weights;
end

% Computing conditional standard deviation, dropping initial window
ConditionalVol      = sqrt(PortCondVarDCC);

% Calculate MSE to compare forecasting accuracy 
for i = 1:DataLength-Window
    Diff             = 0;
    Diff             = Diff+(PortCondVarDCC(i,1)- PortRetSq(i+253,1))^2 ;
end
DCCMSE               = Diff/DataLength ;


%% *** Part 18 - Graphing All Of The Results ***
%  *** in-sample univariate GARCH and rolling correlations ***                            
% Produce charts comparing output 
% Creating time labels
StartDate             = 2006;
EndDate               = 2009;
TimeLabels            = linspace(StartDate, EndDate, DataLength - Window);

% Produce figure
figure(a);
% MSFT
titleStr = 'Volatility Estimates for MSFT';
title(titleStr);
xlabel('Dates');
axis([StartDate EndDate 0 1]);
grid on;
hold on;
plot(TimeLabels, sqrt(252) * MyDataSqr((Window+1:end),1), 'LineStyle', '-' , ...
                  'LineWidth', 1, 'Color', 'red');
plot(TimeLabels, sqrt(252) * MSFTGARCHVol, 'LineStyle', '-' , ...
                  'LineWidth', 1, 'Color', 'magenta');
plot(TimeLabels, sqrt(252) * MSFTGJRGARCHVol, 'LineStyle', '-' , ...
                  'LineWidth', 1, 'Color', 'blue');
plot(TimeLabels, sqrt(252) * MSFTEGARCHVol, 'LineStyle', '-' , ...
                  'LineWidth', 1, 'Color', 'green');
plot(TimeLabels, sqrt(252) * CCCVolMSFT, 'LineStyle', '--', ...
                 'LineWidth', 1, 'Color', 'yellow');
plot(TimeLabels, sqrt(252) * DCCVolMSFT, 'LineStyle', '--', ...
                 'LineWidth', 1, 'Color', 'black');
legend('Actual Data', 'GARCH', 'GJR-GARCH','EGARCH','CCC','DCC','location','best');              
set(gcf, 'color', 'white'); 
ylabel('Volatility');

a=a+1;


 % XOM 
% Produce figure
figure(a);
titleStr = 'Volatility Estimates for EXXON';
title(titleStr);
xlabel('Dates');
axis([StartDate EndDate 0 1.5]);
grid on;
hold on;
plot(TimeLabels, sqrt(252) * MyDataSqr((Window+1:end),2), 'LineStyle', '-' , ...
                  'LineWidth', 1, 'Color', 'red');
plot(TimeLabels, sqrt(252) * XOMGARCHVol, 'LineStyle', '-' , ...
                  'LineWidth', 1, 'Color', 'magenta');
plot(TimeLabels, sqrt(252) * XOMGJRGARCHVol, 'LineStyle', '-' , ...
                  'LineWidth', 1, 'Color', 'blue');
plot(TimeLabels, sqrt(252) * XOMEGARCHVol, 'LineStyle', '-' , ...
                  'LineWidth', 1, 'Color', 'green');
plot(TimeLabels, sqrt(252) * CCCVolXOM, 'LineStyle', '--', ...
                 'LineWidth', 1, 'Color', 'yellow');
plot(TimeLabels, sqrt(252) * DCCVolXOM, 'LineStyle', '--', ...
                 'LineWidth', 1, 'Color', 'black');
legend('Actual Data', 'GARCH', 'GJR-GARCH','EGARCH','CCC','DCC','location','best');              
set(gcf, 'color', 'white'); 
ylabel('Volatility');

a=a+1;

 % C
% Produce figure
figure(a);
titleStr = 'Volatility Estimates for C';
title(titleStr);
xlabel('Dates');
axis([StartDate EndDate 0 1.5]);
grid on;
hold on;
plot(TimeLabels, sqrt(252) * MyDataSqr((Window+1:end),3), 'LineStyle', '-' , ...
                  'LineWidth', 1, 'Color', 'red');
plot(TimeLabels, sqrt(252) * CGARCHVol, 'LineStyle', '-' , ...
                  'LineWidth', 1, 'Color', 'magenta');
plot(TimeLabels, sqrt(252) * CGJRGARCHVol, 'LineStyle', '-' , ...
                  'LineWidth', 1, 'Color', 'blue');
plot(TimeLabels, sqrt(252) * CEGARCHVol, 'LineStyle', '-' , ...
                  'LineWidth', 1, 'Color', 'green');
plot(TimeLabels, sqrt(252) * CCCVolC, 'LineStyle', '--', ...
                 'LineWidth', 1, 'Color', 'yellow');
plot(TimeLabels, sqrt(252) * DCCVolC, 'LineStyle', '--', ...
                 'LineWidth', 1, 'Color', 'black');
legend('Actual Data', 'GARCH', 'GJR-GARCH','EGARCH','CCC','DCC','location','best');                
set(gcf, 'color', 'white'); 
ylabel('Volatility');
a=a+1;


% Correlation Plot
% Produce figure
figure(4)
% MSFT XOM 
subplot(3, 1, 1);
titleStr = 'Correlation Estimates for  MSFT and XOM  ';
title(titleStr);
xlabel('Dates');
axis([StartDate EndDate -1.2 1.2]);
grid on;
hold on;
plot(TimeLabels, MultCorr(:,1), 'LineStyle', '-' , ...
                  'LineWidth', 1, 'Color', 'blue');
plot(TimeLabels, DCCCorr12, 'LineStyle', '--', ...
                 'LineWidth', 1, 'Color', 'red');
legend('20-Day', 'Out-of-Sample DCC', 'location','best');              
set(gcf, 'color', 'white'); 
ylabel('Correlation');

% MSFT C
subplot(3, 1, 2);
titleStr = 'Correlation Estimates for  MSFT and C  ';
title(titleStr);
xlabel('Dates');
axis([StartDate EndDate -1.2 1.2]);
grid on;
hold on;
plot(TimeLabels, MultCorr(:,2), 'LineStyle', '-' , ...
                  'LineWidth', 1, 'Color', 'blue');
plot(TimeLabels, DCCCorr13, 'LineStyle', '--', ...
                 'LineWidth', 1, 'Color', 'red');
legend('20-Day', 'Out-of-Sample DCC', 'location','best');              
set(gcf, 'color', 'white'); 
ylabel('Correlation');

% XOM C
subplot(3, 1, 3);
titleStr = 'Correlation Estimates for  XOM and C  ';
title(titleStr);
xlabel('Dates');
axis([StartDate EndDate -1.2 1.2]);
grid on;
hold on;
plot(TimeLabels, MultCorr(:,3), 'LineStyle', '-' , ...
                  'LineWidth', 1, 'Color', 'blue');
plot(TimeLabels, DCCCorr23, 'LineStyle', '--', ...
                 'LineWidth', 1, 'Color', 'red');
legend('20-Day', 'Out-of-Sample DCC', 'location','best');              
set(gcf, 'color', 'white'); 
ylabel('Correlation');

a=a+1;

%% *** Part 19 - Performing VaR analysis with portfolio ***
% Parameters for VaR backtest
p                   = 0.01;            % VaR level
InvCDF              = norminv(p);
ValueAtRiskEst      = - InvCDF * ConditionalVol;

% Finding breaches
Breach              = (PortLogRet((Window+1):end, :) < - ValueAtRiskEst);
LenBreach           = size(Breach, 1);
pHat                = sum(Breach)/size(Breach, 1);
VRtest              = pHat/p;


% Unconditional coverage test
V1                  = sum(Breach);
V0                  = size(Breach, 1) - V1;

LogLConst           = V1 * log(p) + V0 * log(1-p);
LogLUnconst         = V1 * log(pHat) + V0 * log(1-pHat);
TestStatisticUC     = -2 * (LogLConst - LogLUnconst);
SignificanceUC      = 1 - chi2cdf(TestStatisticUC, 1);

% Conditional coverage test
% Estimating the transition probabilities
Aux01               = (ones(LenBreach - 1, 1) - Breach(1:(end - 1), :)) ...
                         .* Breach(2:end, :);
p01                 = sum(Aux01) / (LenBreach - 1 ...
                         - sum(Breach(1:(end -1), :)));
                     
Aux11               = Breach(1:(end - 1), :) .* Breach(2:end, :);
p11                 = sum(Aux11) / sum(Breach(1:(end -1), :));

% Computing log-likelihood; note that the constrained version here is 
% equal to the unconstrained version above
LogLConst           = LogLUnconst;

% Running loop to compute the unconstrained log-likelihood
LogLUnconst         = Breach(1, 1) * log(pHat) + ...
                         (1 - Breach(1, 1)) * log(1 - pHat);
for i=2:size(Breach, 1)
    if     (Breach(i-1, 1) == 0 && Breach(i, 1) == 0)
        LogLUnconst = LogLUnconst + log(1 - p01);
    elseif (Breach(i-1, 1) == 0 && Breach(i, 1) == 1)
        LogLUnconst = LogLUnconst + log(p01);
    elseif (Breach(i-1, 1) == 1 && Breach(i, 1) == 0)
        LogLUnconst = LogLUnconst + log(1 - p11);
    elseif (Breach(i-1, 1) == 1 && Breach(i, 1) == 1)
        LogLUnconst = LogLUnconst + log(p11);
    end
end

TestStatisticCC     = -2 * (LogLConst - LogLUnconst);
SignificanceCC      = 1 - chi2cdf(TestStatisticCC, 1);

% Calculate VaR at 1%,5% and 10%
alpha               = [0.01,0.05,0.1]; % VaR level
InvCDF              = norminv(alpha);
VaRDCC              = NaN(size(ConditionalVol,1),3);
for i = 1:3
VaRDCC(:,i)         = - InvCDF(1,i) * ConditionalVol;
end










