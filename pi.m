%% Find the value of pi using simulations

clear all;
clc;
close all;

N=100;
S=rand(N,2);
C=((S(:,1).^2+S(:,2).^2)<=1);
PI=4*sum(C)/N;
% disp(PI);    
    
plot(S(:,1),S(:,2),'b.')    
hold on
plot(C.*S(:,1),C.*S(:,2),'r.')
hold on
plot([0:0.0001:1],sqrt(1-[0:0.0001:1].^2),'k')