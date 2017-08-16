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
