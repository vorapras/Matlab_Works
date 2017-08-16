
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
                  'LineWidth', 1, 'Color', 'cyan');
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
                  'LineWidth', 1, 'Color', 'cyan');
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
                  'LineWidth', 1, 'Color', 'cyan');
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


