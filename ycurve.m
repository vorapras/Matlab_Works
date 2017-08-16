function x=ycurve(yieldgrid,yielddate,calendar,startdate,CF,price,N,T);

%first interpolate to get yield for every date in dataset
for t=1:T;
    [a b]=min(abs(yielddate-calendar(t))+10000*(yielddate>=calendar(t)));
    if b==length(yielddate); b=b-1; end;
    mult=(calendar(t)-yielddate(b))/(yielddate(b+1)-yielddate(b));
    yield(t,1)=yieldgrid(b)+mult*(yieldgrid(b+1)-yieldgrid(b));
end;
%next calculate the discounted present value of CF and compare to price
for i=1:N;
    err(i)=price(i)-sum(CF(:,i)./(yield.^(calendar([1:T]')-startdate)));
end;

x=mean(abs(err));