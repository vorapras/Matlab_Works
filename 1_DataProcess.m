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
Data_D           = importdata('Thanachai_Data.csv');
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

