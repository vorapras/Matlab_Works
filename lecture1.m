%% Code for FM457A and FM457B - Lecture 1

%always start your code by clearing the memory
clear all;
close all;
clc;

%define a variable x, with value 1
x=1;

%define a variable y, with value 5
y=5;

%define a 1x5 vector V1
V1=[1.1 .99 1.07 .87 .95];

%define a 1x5 vector V2
V2=[.9 1.03 1.03 1.05 .99];

%see which variables exist, use:
whos;

%delete variables x and y 
clear x y;

%define a 3x5 matrix M1
M1=[1.01 1.01 1.03 .99 .93
    .99 .98 .97 .96 .95
    1.21 .77 .89 1.13 1.05];

%find out size of vectors or matrices
size(M1)

%Define a 3x5 matrix M2, which uses V1 and V2 as pieces.
%Note that I can use semicolon to start the next line 
%instead of actually starting on the next line
M2=[V1; V2; V1];

%Define a 3x5 matrix M3 with all of its entries being zero
M3=zeros(3,5);

%Define a 3x5 matrix M4 with all of its entries being ones
M4=ones(3,5);

%Define a variable z, whose value is the 3rd row, 2nd column of M1
z=M1(3,2);

%% Elementary Algebra

%violating rules of matrix algebra
M5 = M1 + V1;

%matrix addition
M5 = M1 + M2

%element by element multiplication
M5 = M1.*M2

%element by element division
M5 = M1./M2

%transpose a matrix
M5 = M1'

%matrix multiplication
M6 = M5*M1

clear M5

%scalar mutliplication
x=2;
M5 = M1*x


%% Ranges of numbers

%number ranges
V3 = [1:5]

%number ranges with different step size
V3 = [1:3:13]

M5 = M1(1:2,1:5)

V3 = M1(:,2)

%% Descriptive Statistics - Useful Functions

%load data
data_bp;
data_vlcm;

%look at part of the data
bp(1:5,:)

%maximum of realized BP returns
max(bp(:,4))

%minimum of realized Volcom returns and location of minimum
[a b] = min(vlcm(:,4))

%mean of daily returns in percent
100*mean([bp(:,4) vlcm(:,4)])

%std of daily returns in percent
100*std([bp(:,4) vlcm(:,4)])

%correlation coefficient of return time series
corrcoef([bp(:,4) vlcm(:,4)])

%% Plot Time Series

%capture number of rows and columns of a vector
[c d] = size(vlcm);

%create a time index
time = 1:c;

%plot Volcom's return time series
plot(time,vlcm(:,4),'b');

%tell Matlab to hold on to existing plot
hold on;

%plot BP's return time series in the same graph
plot(time,bp(:,4),'r--');

hold off; %turns off hold on

%multiple plots
subplot(2,1,1); %splits plot into 2 plots
%split is horizontal, currently we plot on the top

%plot Volcom's price time series in green
plot(time, vlcm(:,3),'g');

%label the axes of the first plot
title('Volcom Price');
xlabel('Time');
ylabel('Return');

subplot(2,1,2); %now we plot below

%scatterplot of Volcom's and BP's return time series
plot(vlcm(:,4),bp(:,4),'k.');

%label the axes of the second plot
title('Scatter Plot of Returns');
xlabel('Volcom');
ylabel('British Petroleum');

%% Saving the work

%saves the variables x and z in a file called workspace1.mat
save workspace1 x z

%saves all variables in a file called workspace1.mat
save workspace1

%clear workspace
clear all;
close all;
clc;

%loads the variables in workspace1.mat
load workspace1 























