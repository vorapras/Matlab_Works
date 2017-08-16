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