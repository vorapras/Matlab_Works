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
PortCondVol         = PortRetEGARCHVol ; %PortRetGARCHVol ; % PortRetGJRGARCHVol ; % PortRetEGARCHVol ;
p                   = 0.1;            % VaR level
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
