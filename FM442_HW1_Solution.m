%% London School of Economics
%  FM442 - Michaelmas Term 2016
%  Author: Domingos Romualdo
%  Date: 11 October 2016
%  Assignment 1 - Suggested Solution

%% *** Clearing workspace *** 
clc;
clear all;
clf;

%% *** Part 1: Loading monthly data ***
Data_M = importdata('FM442 - HW1 - Monthly Data.csv');

% Note: since the original file contains text and data, the import will
% have two parts:
%
% - Data_M.data contains numbers
% - Data_M.textdata contains text (the first row containing all the labels)
%
% Note that Data_M.textdata has one fewer row than Data_M.data 

%% *** Part 2 - Cleaning up and organizing the data ***
% Retrieving information about data
[NRows, NCols] = size(Data_M.data);

% Preallocating space for matrices with PermNo and Dates
PermNo = NaN(NRows, 1);
Dates  = NaN(NRows, 1);

for i = 1:NRows
    % Note: the str2num function converts text to data
    PermNo(i,1) = str2num(Data_M.textdata{i+1, 1});
    Dates( i,1) = str2num(Data_M.textdata{i+1, 2});
end

% Collecting data in a single matrix
DataNum       = [Dates, PermNo, Data_M.data];

% Getting unique values of stock identifiers and dates
PermNoUnique  = unique(PermNo);
DatesUnique   = unique(Dates);

% Computing number of stocks and dates
NStocks       = length(PermNoUnique);
NDates        = length(DatesUnique);

%% *** Part 3 - Replicating returns data ***

% First, we try to replicate the price returns (RETX) using the original
% prices (PRC) and the adjustment factor (CFCPR).  Note that original
% prices are in column 5, adjustment factors in column 9, and returns (ex-
% dividends) in column 10 of DataNum.  We will do so on a 
% security-by-security basis.

% Allocate memory for vector that we will use to check whether we can
% replicate returns
MaxDiff          = NaN(NStocks, 1);

for i = 1:NStocks
   % Extract the rows that correspond to that stock
   Indices       = find(PermNo==PermNoUnique(i));
   StockData     = DataNum(Indices, :);
   
   % Sort the data according to date
   StockData     = sortrows(StockData);
   
   % Find out which rows contain dates that are not equal to the one in
   % the next row
   Dates1        = StockData(1:(end-1), 1);
   Dates2        = StockData(2:end, 1);
   Indices2      = find(Dates2 ~= Dates1);
   
   % Check if date in last row is unique, and if so include it in sample
   if (StockData(end - 1, 1) ~= StockData(end, 1))
       Indices2  = [Indices2; size(StockData, 1)];
   end
   
   % Extract data for this stock for those dates only
   StockData     = StockData(Indices2, :);
   
   % Perform computation to replicate returns; note that original prices
   % are in column 5, adjustment factors in column 9, and returns (ex-
   % dividends) in column 10 of DataNum
   AdjPrice      = StockData(:, 5) ./ StockData(:, 9);
   ComputedRet   = AdjPrice(2:end, :) ./ AdjPrice(1:(end-1), :) - 1;
   AbsDiff       = abs(ComputedRet - StockData(2:end, 10));
   
   % Compute maximum value of the difference, store it in MaxDiff
   MaxDiff(i, 1) = max(AbsDiff);
end

% Next, we try to replicate the returns (RET) using the price returns
% prices (RETX), dividend amounts (DIVAMT) and the factors to adjust
% dividends (FACPR).  Note that these are in columns 6, 10, 3 and 4 of 
% DataNum.  Since we have already replicated the ex-dividend returns with
% prices, we will use these instead of PRC.

% There are a number of steps we need to go through here.  First, since
% dividends are sometimes reported in multiple rows, we must consolidate
% them for each stock and each date.  This will involve running a loop,
% which, as it turns out, works better when executed in reverse order.
% When that is completed, we will compute the returns with dividends using
% adjustment factors.

% Note: in the process of replicating this data, we will also create
% matrices that contain total returns for each stock, as well as the S&P
% 500 and the value-weighted return for this universe

% Allocate memory for vector that we will use to check whether we can
% replicate returns
MaxDiff2         = NaN(NStocks, 1);

% Allocate memory for total returns for each stock
TR_Stocks        = NaN(NDates, NStocks);

for i = 1:NStocks
   % Extract the rows that correspond to that stock
   Indices       = find(PermNo==PermNoUnique(i));
   StockData     = DataNum(Indices, :);
   
   % Sort the data according to date
   StockData     = sortrows(StockData);
   
   % Computing number of rows in our data
   NRows         = size(StockData, 1);
   
   % Running loop to consolidate data into single rows per date
   for j=NRows:-1:2
      if(StockData(j, 1) == StockData(j - 1, 1))
         StockData(j - 1, 3) = StockData(j - 1, 3) + StockData(j, 3);
         StockData(j, :) = [];
      end 
   end
   
   % Now, compute total returns from data
   % One minor difficulty is that many dividend numbers are missing,
   % so Matlab treats them as NaN.  In order to add those, we need to
   % replace them by zeros, which is done in the next two lines
   AuxDiv        = StockData(:, 3);
   AuxDiv(isnan(AuxDiv)) = 0; 
   Numerator     = (AuxDiv + StockData(:, 5)) ./ StockData(:, 9);
   Denominator   = StockData(:, 5) ./ StockData(:, 9);
   ComputedRet   = Numerator(2:end, :) ./ Denominator(1:(end-1), :) - 1;
   AbsDiff       = abs(ComputedRet - StockData(2:end, 6));
   
   % Compute maximum value of the difference, store it in MaxDiff
   MaxDiff2(i, 1) = max(AbsDiff);
   
   % Store total returns for this stock in matrix
   TR_Stocks(:, i) = StockData(:, 6);
   
   % Store total returns for SPX and value-weighted index
   if i == 1
       TR_SPX      = StockData(:, 15);
       TR_VW       = StockData(:, 11);
   end
end

% The MaxDiff will be small for any stock which has only had ordinary
% dividends and stock splits/reverse splits.  In the case of spinoffs,
% additional information is required to compute total returns, so we can't
% reproduce the data.

%% *** Part 4: Working with Monthly Data ***

% Compiling single matrix with all total returns data
TR                 = [TR_Stocks, TR_SPX, TR_VW];  

% Computing total returns indices
TRIndices          = NaN(NDates + 1, NStocks + 2);
TRIndices(1, :)    = ones(1, NStocks + 2);

for i=1:NDates
    TRIndices(i+1, :) = TRIndices(i, :) .* (1 + TR(i, :));
end

% Computing log returns
LogReturns         = log(1 + TR);

% Computing statistics
% Note the conversion of mean and standard deviation into annual units
MonthlyStats       = NaN(4, NStocks + 2);
MonthlyStats(1, :) = 12 * mean(LogReturns);
MonthlyStats(2, :) = sqrt(12) * std(LogReturns);
MonthlyStats(3, :) = skewness(LogReturns);
MonthlyStats(4, :) = kurtosis(LogReturns);

% Producing charts for total returns of each of the stocks versus the
% value-weighted index

% First, create time labels using linspace
StartDate     = 1989;
EndDate       = 2015; 
TimeLabels    = linspace(StartDate, EndDate, NDates + 1);

% Creating list of tickers
% Note that strings must be the same length, or Matlab will be unhappy
Tickers       = ['MSFT'; 'XOM '; 'GE  '; 'JPM '; 'INTC'; 'C   '];

% These contain the appropriate upper limits for the y-axis
AuxVec        = [140; 20; 20; 20; 80; 40];

% In order to create multiple figures, must set figure number
FigNo         = 1;

% Loop to create all the figures
for i = 1:NStocks
    % Parameters for chart
    ChartData   = [TRIndices(:, i) TRIndices(:, 8)];
    Stock_Str   = Tickers(i, :);
    Axis1       = [StartDate EndDate 0 AuxVec(i)];
    
    % Creating figure
    figure(FigNo)
    title(Stock_Str);
    ylabel('Total Returns Index');
    xlabel('Dates');
    axis(Axis1);
    grid on;
    hold on;
    plot(TimeLabels, ChartData(:,1), 'LineStyle', '-' , ...
                      'LineWidth', 1, 'Color', 'black');
    plot(TimeLabels, ChartData(:,2), 'LineStyle', '--', ...
                     'LineWidth', 1, 'Color', 'blue');
    legend1 = legend('Total Return', 'Value-Weighted Returns', ...
                     'location','best');
    set(gcf, 'color', 'white');

    FigNo = FigNo + 1;
end

%% *** Part 5: Working with Daily Returns ***

% Import data from Excel spreadsheet
% Here, I use the xlsread function
Data_D_CRSP       = xlsread('FM442 - HW1 - Daily Data.xlsx', ...
                     'HPR Daily', 'A1:D6554');

Data_D_Yahoo      = xlsread('FM442 - HW1 - Daily Data.xlsx', ...
                     'Yahoo Data - Import', 'A4:C6746');

% Delete the column with identifiers in the CRSP data
Data_D_CRSP(:, 1) = [];

% Sort the data by the dates in both cases
Data_D_CRSP    = sortrows(Data_D_CRSP);
Data_D_Yahoo   = sortrows(Data_D_Yahoo);

% Delete Yahoo data past the end of 2015
% Note that once these lines are executed, the CRSP and Yahoo data have the
% same number of rows
YahooDates     = Data_D_Yahoo(:, 1);
Indices        = find(YahooDates > 20151231);
Data_D_Yahoo(Indices, :) = [];

% The Yahoo data is a total returns index, while the CRSP is a returns
% series.  Thus, construct total returns series for CRSP
NumRows        = size(Data_D_CRSP, 1);
TRI_CRSP       = NaN(NumRows + 1, 2);
TRI_CRSP(1, :) = ones(1, 2);

for i=1:NumRows
    TRI_CRSP(i+1, :) = TRI_CRSP(i, :) .* (1 + Data_D_CRSP(i, 2:3));
end

% Delete first row, normalize data for CRSP
TRI_CRSP(1, :) = [];
TRI_CRSP       = TRI_CRSP ./ repmat(TRI_CRSP(1, :), NumRows, 1);

% Construct and normalize Yahoo data
TRI_Yahoo      = Data_D_Yahoo(:, 2:3);
TRI_Yahoo      = TRI_Yahoo ./ repmat(TRI_Yahoo(1, :), NumRows, 1);

% Produce charts to compare the data
NDates         = size(TRI_Yahoo, 1);
TimeLabels    = linspace(StartDate, EndDate, NDates);

% Creating list of tickers
% Note that strings must be the same length, or Matlab will be unhappy
Tickers       = ['MSFT'; 'SPX '];

% These contain the appropriate upper limits for the y-axis
AuxVec        = [140; 8];

for i = 1:2
    % Parameters for chart
    ChartData   = [TRI_CRSP(:, i) TRI_Yahoo(:, i)];
    Stock_Str   = Tickers(i, :);
    Axis1       = [StartDate EndDate 0 AuxVec(i)];
    
    % Creating figure
    figure(FigNo)
    title(Stock_Str);
    ylabel('Total Returns Index');
    xlabel('Dates');
    axis(Axis1);
    grid on;
    hold on;
    plot(TimeLabels, ChartData(:,1), 'LineStyle', '-' , ...
                      'LineWidth', 1, 'Color', 'black');
    plot(TimeLabels, ChartData(:,2), 'LineStyle', '--', ...
                     'LineWidth', 1, 'Color', 'blue');
    legend1 = legend('CRSP Series', 'Yahoo Series', ...
                     'location','best');
    set(gcf, 'color', 'white');

    FigNo = FigNo + 1;
end

% Compute log returns, statistics
ReturnsDaily       = Data_D_CRSP(:, 2:3);
LogReturnsD        = log(1 + ReturnsDaily);

% Computing statistics
% Note the conversion of mean and standard deviation into annual units
DailyStats         = NaN(4, 2);
DailyStats(1, :)   = 252 * mean(LogReturnsD);
DailyStats(2, :)   = sqrt(252) * std(LogReturnsD);
DailyStats(3, :)   = skewness(LogReturnsD);
DailyStats(4, :)   = kurtosis(LogReturnsD);

% Extracting Monthly Index from the daily data to compare with the
% monthly series.  The key is to identify rows that correspond to end-
% month dates.  First, compute the month for each observation
DailyDates         = Data_D_CRSP(:, 1);
AuxSeries          = (DailyDates - mod(DailyDates, 100))/100;
DateMonth          = mod(AuxSeries, 100);

% Next, determine which observations are month-end
MonthEndRows       = find(DateMonth(1:(end-1), :) ~= DateMonth(2:end, :));
MonthEndRows       = [MonthEndRows; size(DailyDates, 1)];
DailyIndices       = TRI_CRSP(MonthEndRows, :);

% Normalize month-end index by first row
DailyIndices       = DailyIndices ./ repmat(DailyIndices(1, :), ...
                        size(DailyIndices, 1), 1);
                    
% Get info from monthly returns for MSFT and SPX (columns 1 and 7)
% Since monthly data started in Dec 1989, need to delete first two rows
% Also, renormalize series
ColumnIndices      = [1, 7];
MonthlyIndices     = TRIndices(3:end , ColumnIndices);
MonthlyIndices     = MonthlyIndices ./ repmat(MonthlyIndices(1, :), ...
                        size(MonthlyIndices, 1), 1);

% Set options for charts to compare the data
NDates             = size(DailyIndices, 1);
TimeLabels         = linspace(StartDate, EndDate, NDates);

% Produce charts
for i = 1:2
    % Parameters for chart
    ChartData   = [DailyIndices(:, i) MonthlyIndices(:, i)];
    Stock_Str   = Tickers(i, :);
    Axis1       = [StartDate EndDate 0 AuxVec(i)];
    
    % Creating figure
    figure(FigNo)
    title(Stock_Str);
    ylabel('Total Returns Index');
    xlabel('Dates');
    axis(Axis1);
    grid on;
    hold on;
    plot(TimeLabels, ChartData(:,1), 'LineStyle', '-' , ...
                      'LineWidth', 1, 'Color', 'black');
    plot(TimeLabels, ChartData(:,2), 'LineStyle', '--', ...
                     'LineWidth', 1, 'Color', 'blue');
    legend1 = legend('Daily Series', 'Monthly Series', ...
                     'location','best');
    set(gcf, 'color', 'white');

    FigNo = FigNo + 1;
end
