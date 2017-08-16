%% London School of Economics
%  FM442 - Michaelmas Term 2016
%  Author: Domingos Romualdo
%  Date: 15 November 2016
%  Assignment 5 - Suggested Solution

% *** Preliminaries ***
clc;
clear all;
clf;
close all;
FigNo            = 1;

%% *** Part 2: Random Number generation and analysis ***

%% Item a: Uniformly distributed numbers
Size             = 100000;
SampleUnif       = rand(Size, 1);

% Computing properties of the data
UnifMean         = mean(SampleUnif);   % Expect close to 0.5
UnifStdDev       = std(SampleUnif);    % Expect close to 1/sqrt(12) = 0.288

% Producing histogram of the data (analog to pdf)
figure(FigNo);
histogram(SampleUnif);
TitleStr         = 'Histogram of Uniformly Distributed Data';
title(TitleStr);
set(gca, 'XTickMode', 'manual');  
FigNo            = FigNo + 1;

% Charting the empirical cumulative distribution function
cdfplot(SampleUnif);
titleStr = 'Cumulative Distribution of Uniformly Distributed Data';
title(titleStr);
FigNo            = FigNo + 1;

%% Item b: Normally distributed random numbers with randn
SampleNormal     = randn(Size, 1);

% Computing properties of the data
NormalMean       = mean(SampleNormal);   % Expect close to 0
NormalStdDev     = std(SampleNormal);    % Expect close to 1

% Producing histogram of the data (analog to pdf)
figure(FigNo);
histogram(SampleNormal);
TitleStr         = 'Histogram of Normally Distributed Data';
title(TitleStr);
set(gca, 'XTickMode', 'manual');  
FigNo            = FigNo + 1;

% Charting the empirical cumulative distribution function
cdfplot(SampleNormal);
titleStr = 'Cumulative Distribution of Normally Distributed Data';
title(titleStr);
FigNo            = FigNo + 1;

%% Item c: Normally distributed random numbers with norminv
SampleNormal2     = norminv(rand(Size, 1));

% Computing properties of the data
NormalMean2      = mean(SampleNormal2);   % Expect close to 0
NormalStdDev2    = std(SampleNormal2);    % Expect close to 1

% Producing histogram of the data (analog to pdf)
figure(FigNo);
histogram(SampleNormal2);
TitleStr         = 'Histogram of Normally Distributed Data';
title(TitleStr);
set(gca, 'XTickMode', 'manual');  
FigNo            = FigNo + 1;

% Charting the empirical cumulative distribution function
cdfplot(SampleNormal2);
titleStr = 'Cumulative Distribution of Normally Distributed Data';
title(titleStr);
FigNo            = FigNo + 1;

%% Item d: t(4)distributed random numbers
SampleTDist      = tinv(rand(Size, 1), 4);

% Computing properties of the data
TDistMean        = mean(SampleTDist);   % Expect close to 0
TDistStdDev      = std(SampleTDist);    % Expect close to 1.414

% Producing histogram of the data (analog to pdf)
figure(FigNo);
histogram(SampleTDist);
TitleStr         = 'Histogram of t(4)-Distributed Data';
title(TitleStr);
set(gca, 'XTickMode', 'manual');  
FigNo            = FigNo + 1;

% Charting the empirical cumulative distribution function
cdfplot(SampleTDist);
titleStr = 'Cumulative Distribution of t(4)-Distributed Data';
title(titleStr);
FigNo            = FigNo + 1;

%% Item e: Reproducible random numbers
% Generate some random data
rng(1);
SampleUnif1      = rand(Size, 1);

% Do some work, generate new sample, check that they are different
Aux1             = mean(SampleUnif1);
Aux2             = std(SampleUnif1);
SampleUnif2      = rand(Size, 1);
TestVal          = max(abs(SampleUnif1 - SampleUnif2));

% Reset the seed, generate new random numbers,
% then verify that they are the same as SampleUnif1
rng(1);
SampleUnif3      = rand(Size, 1);
TestVal          = max(abs(SampleUnif1 - SampleUnif3));

%% *** Part 3: Black-Scholes Options Prices ***

% Setting parameters
P                = 40;         % Current stock price
K                = 40;          % Options strike price
r                = 0.05;        % Interest rate per annum
T                = 1;           % Time to expiration in years
Sigma            = 0.20;        % Stock price volatility

%% Item f: Computing options prices using formula
[CallFormula, PutFormula] = blsprice(P, K, r, T, Sigma);

%% Item g: conducting Monte Carlo simulation to find prices
% Generating random numbers with uniform normal distribution
Size             = 100000;
Eps              = randn(Size, 1);

% Computing stock price at expiration for each realization of Eps
Exponent         = (r - Sigma^2/2) * T * ones(size(Eps)) + ... 
                      Sigma * sqrt(T) * Eps;
StockVal         = P * exp(Exponent);

% Computing call values at expiration
CallValues = StockVal - K*ones(Size, 1);
CallValues(CallValues < 0) = 0;

% Computing put values at expiration
PutValues  = K*ones(Size, 1) - StockVal;
PutValues(PutValues < 0) = 0;

% Taking expected values, discounting them back to the present
CallMonteCarlo   = exp(-r*T) * mean(CallValues);
PutMonteCarlo    = exp(-r*T) * mean(PutValues);
