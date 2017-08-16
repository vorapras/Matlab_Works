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


