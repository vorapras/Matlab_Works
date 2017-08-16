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







