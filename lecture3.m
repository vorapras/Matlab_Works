%% Code for FM457A and FM457B - Lecture 3

% clear workspace
clear all;
close all;
clc;


%% Histograms
X=[2*ones(3,1); 3*ones(5,1); 7*ones(4,1)];
subplot(2,1,1); 
hist(X); %draws histogram of X
subplot(2,1,2); 
hist(X,[0:.25:10]); %draws histogram of X
%on the interval [0 10] with bins of size 0.25
         
%draw random numbers with rand()                  
T=1000;
X=rand(T,1); %creates a matrix of size
%(Tx1) of uniform rv’s on (0,1)
a=5; b=50;
Y=a+b*rand(T,1); %creates a matrix of
%size (Tx1) of unifrom rv’s on (a,a+b)  
subplot(3,1,1);
hist(X,[-.25:.1:1.25]); %draws histogram of X
subplot(3,1,2); 
hist(Y,[.9*a:(b/50):1.1*(a+b)]); %draw
T=100000; Z=a+b*rand(T,1); 
subplot(3,1,3);
hist(Z,[.9*a:(b/50):1.1*(a+b)]); subplot;

%% Law of Large Numbers

%visualize a law of large numbers
X=rand(5,1); disp(mean(X));
X=rand(10,1); disp(mean(X));
clear Y;
for i=1:200;
    Y(i)=mean(rand(i,1));
end;
plot(Y); 

%create discrete random variables
%x=1 with probability p and 0 with (1-p)
p=.25; 
if rand(1,1)<p; 
    x=1; 
else; 
    x=0; 
end;

%x=3 with p=.25, 2 with p=.5, 1 with p=.25
y=rand(20,1); 
x=3.*(y>.75)+2*(y<=.75 & y>.25)+1*(y<=.25);
hist(x,[0:.25:4]);


%% Central Limit Theorem

A=[1 5 10 25 50 100]; %number of draws
for k=1:6;
    j=A(k);   
    for i=1:5000; 
        x=rand(j,1);
        y=(x<.3)*0+(x>=.3)*1; %think of y as a biased coin flip
        z(i,k)=mean(y);
    end;
end;

%plot the different scenarios
for k=1:6;
    subplot(3,2,k); 
    hist(z(:,k),50);
end;
subplot;

%draw normal random variables with randn()
X=randn(10000,1);
subplot(2,1,1); hist(X,50);

%in order to adjust mean and/or variance, just skew and shift
m=1.1; s=.16; X=m+s*randn(10000,1);
subplot(2,1,2); hist(X,50);

%a simple securtiy process
T=10000; mu=1.1; sigma=.3; 
x=randn(T,1);
R=mu+sigma*x;
subplot; hist(R,50);

%Question: Are there any problems with this?
%Answer: LL


%% Calibration

%calculate calibration parameters Microsoft
data_msft;
disp(mean(msft(:,4)));
disp(std(msft(:,4)));
%Microsoft returns have a daily mean of
%.097% and standard deviation of 2.21%
subplot(2,1,1); hist(msft(:,4),[-.2:.01:.2]);
axis([-.2 .2 0 800]); xlabel('Actual');

%simulate Microsoft
T=3022; 
mu=.00097-.5*.0221^2; sigma=.0221;
x=randn(T,1);
r=exp(mu+sigma*x)-1; 
subplot(2,1,2);
hist(r,[-.2:.01:.2]); axis([-.2 .2 0 800]);
axis([-.2 .2 0 800]); xlabel('Simulated');
disp(mean(r)); disp(std(r));
disp([skewness(r) skewness(msft(:,4))]);
disp([kurtosis(r) kurtosis(msft(:,4))]);

