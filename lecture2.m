%% Code for FM457A and FM457B - Lecture 2

% clear workspace
clear all;
close all;
clc;

%% Logical Operators

%logical statements
1 == 1

1 == 2

x=5; %assigns value 5 to x

x==5; %checks if x is equal to 5, returns either 1 or 0

x=(x==5); %checks if x is equal to 5, then assigns True (1) or False (0) to x


%some examples using logical operators

A=[zeros(3,1); ones(3,1); 2*ones(3,1); 3*ones(3,1)];
A'

in1=(A>0);
in1'
A(in1)'

A(A>0); % is same as A(in1)
A(A>0)'

in2=(A<1 | A>2);
in2'
A(in2)
mean(A(in2))

A'
in3=(A>1 & A<3);
in3'
A(in3)'

A'
in4=(A~=2);
in4'
A(in4)'

%% Control Structures - if statements

%if statement
if A(1)==0;
  x=5; y=x;
end;

%elseif statement
if A(3)==A(4);
    x=A(5); y=A(4);
elseif A(3)==0;
    x=5; y=0;
elseif A(3)==1;
    x=4; y=5;
else;
    x=3; y=8;
end;

%% Control Structures - for and while loops

%for loop
T=100; s=0; x=0;
for i=1:T;
    s=s+i;
    x=x+i*i;
end;

%nested loop
for i=1:5;
    for j=1:5;
        B(i,j)=min(i,j);
    end;
end;

%loop without using for
i=0;
while i<10;
    i=i+1;
    disp(i);
end;

%using matrix operations instead of loops
x=0; 
for i=1:5; 
    x=x+i*i; 
end;
%Alternative:
A=[1:5]; y=sum(A.*A);

%corresponding results
x
y

%load data
data_bp;

%slightly longer example
[T L]=size(bp);
s=0; s2=0;
for i=1:T;
    s=s+bp(i,4);
    s2=s2+(bp(i,4)^2);
end;
M=s/T; StD=sqrt((s2/T)-M*M); 
disp([mean(bp(:,4)) M]);
disp([std(bp(:,4)) StD]);


%% Functions

%apply myfunction function
q = myfunction(3,4)

%apply winsorize function and plot result
X=winsorise(bp(:,4),.05,.95);
plot(bp(:,4));
hold on;
grid on;
plot(X(:,1),'r');

%apply npv function and plot result
d=1; r=.05; g=.02; T=20; Ts=10;
p=npv(d,r,g,T,Ts);
T=200; p=zeros(T,1);
for t=1:T;
    p(t)=npv(d,r,g,t,t);
end;   
hold off;
plot([1:T],p,'b');
hold on;
plot([1:T],ones(1,T)*(d/(r-g)),'r--');

