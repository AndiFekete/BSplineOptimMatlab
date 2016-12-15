function wp=predict_mse(f,g,c,k,t,pp)
x=1:length(f);
m=length(t)-2;
p=pp-1;
if p<1 || p>m
    error('The first and the last knot must be kept!');
end
N=length(f);
tt=[t(1).*ones(1,k-1) t t(end).*ones(1,k-1)];
rrho=[tt(1:k-1+p) tt(k+1+p:end)];
%alpha=(t(p)-rrho(1:end-k))./(rrho(k+1:end)-rrho(1:end-k));
%alpha=[0 alpha 0];
%(t(pp)-rrho(p+1:p+k-1))
%(rrho(p+k:p+2*(k-1))-rrho(p+1:p+k-1))
alpha=(t(pp)-rrho(p+1:p+k-1))./(rrho(p+k:p+2*(k-1))-rrho(p+1:p+k-1));

d=zeros(1,m+k-1);
%d(1:p-1)=c(1:p-1);
d(1:p)=c(1:p);
%d(p-1+k-1:end)=c(p+k-1:end);
d(p-1+k:end)=c(p+k:end);
for i=p+1:p+k-2
%    d(i)=(c(i)-(1-alpha(i))*d(i-1))/alpha(i);
    d(i)=(c(i)-(1-alpha(i-p))*d(i-1))/alpha(i-p);
end
%dd=alpha(p-1+k-1)*d(p-1+k-1)+(1+alpha(p-1+k-1)).*d(p-2+k-1);

% dd=zeros(size(c));
% dd(1:p)=d(1:p);
% for i=p+1:p+k-1
%     dd(i)=alpha(i-p)*d(i)+(1-alpha(i-p)).*d(i-1);
% end
% dd(1:p+k-1)
% c(1:p+k-1)



%dd=alpha(end)*d(p-1+k)+(1-alpha(end)).*d(p-1+k-1);

dd=alpha(end)*d(p-1+k)+(1-alpha(end)).*d(p-1+k-1);


dzeta1=(f-g);
%dzeta1=dzeta1+(c(p-1+k-1)-dd)*Bspline(k-1,p-1+k-1,tt,x);
dzeta1=dzeta1+(c(p-1+k)-dd)*Bspline(k-1,p-1+k,tt,x);
dzeta1=sum(dzeta1.^2)/N;



d=zeros(1,m+k-1);
d(1:p)=c(1:p);
d(p-1+k:end)=c(p+k:end);
for i=p+k-1:-1:p+2
    d(i-1)=(c(i)-alpha(i-p)*d(i)) / (1-alpha(i-p));
end
%dd=alpha(p)*d(p)+(1-alpha(p)).*d(p-1);
dd=alpha(1)*d(p+1)+(1-alpha(1)).*d(p);

dzeta2=(f-g);
%dzeta2=dzeta2+(c(p)-dd)*Bspline(k-1,p,tt,x);
dzeta2=dzeta2+(c(p+1)-dd)*Bspline(k-1,p+1,tt,x);
dzeta2=sum(dzeta2.^2)/N;

wp=min(dzeta1,dzeta2);
