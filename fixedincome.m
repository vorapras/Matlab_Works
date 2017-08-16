%%
clear all;
clc;
close all;

%load data
data_bonds;
N=length(maturity); %number of bonds
calendar=[2008+(9/12):(1/12):2030]';
T=length(calendar);
startdate=2008+(10/12)+((3/31)/12);

%
CF=zeros(T,N);
%make matrix of cashflows
for i=1:N;
    [a Tmat]=min(abs(maturity(i)-calendar));
    CF(Tmat,i)=100+coupon(i)*.5;
    for t=1:Tmat-1;
        if mod(t,6)==0 & calendar(Tmat-t)>startdate;
            CF(Tmat-t,i)=coupon(i)*.5;
        end;
    end;
end;

% calculate yields
for i=1:N;
    prerrbest=100000000;
    for j=1:1000;
        y=0+.2*(j-1)/(1000-1);
        py=sum(CF(:,i)./((1+y).^(calendar-startdate)));
        prerr=abs(price(i)-py);
        if prerr<prerrbest;
            prerrbest=prerr;
            yieldbest(i)=y;
        end;
    end;
end;
plot(maturity,yieldbest,'.','MarkerSize',5);
xlabel('Maturity'); ylabel('Yield');
   

%motnhs in which we will solve for yield
yieldgriddate=[6 12 18 24 30 36 42 48 60 72 96 120 180 T]';
%actual dates of those months
yielddate=calendar(yieldgriddate);
%initial guess (these are yiels of bonds at various maturities taken from printout)
y0=[1.0123 1.0157 1.0150 1.0158 1.0175 1.0204 1.0228 1.0238 1.027 1.028 1.033 1.037 1.042 1.042]';


y1=fminsearch(@ycurve,y0,[],yielddate,calendar,startdate,CF,price,N,T);
y2=fminsearch(@ycurve,y1,[],yielddate,calendar,startdate,CF,price,N,T);
y3=fminsearch(@ycurve,y2,[],yielddate,calendar,startdate,CF,price,N,T);
%check minimum values:
disp([ycurve(y0,yielddate,calendar,startdate,CF,price,N,T) ...
    ycurve(y1,yielddate,calendar,startdate,CF,price,N,T) ...
    ycurve(y2,yielddate,calendar,startdate,CF,price,N,T) ...
    ycurve(y3,yielddate,calendar,startdate,CF,price,N,T)]);
%plot the yield curve
plot(yielddate,y3,'LineWidth',3);



