function [c,mse,s,Rll,zll]=bspline_coeffs(sig,knots,order,p,tt,Rl,zl,show)
x=1:1:length(sig);
f=preserve(x,sig);
y=preserve(x,f)';
if nargin<4
    E=calcmat(order,x,knots);
    [Qll,Rll]=qr(E);
    zll=Qll'*y(2:end-1);
    c=Rll\zll;
    show=0;
else
    B=calcmat_opt(order,x,tt,p);  
    [c,Rll,zll]=solve_opt(Rl,zl,B,y(2:end-1),order,p);
end

s=zeros(1,length(sig));
tt=[knots(1).*ones(1,order-1) knots knots(end).*ones(1,order-1)];
c=[0; c; 0];
for i=1:1:length(c)
    p=c(i)*Bspline(order-1,i,tt,x);
    s=s+p;
%     plot(x,c(i)*Bspline(order-1,i,tt,x))
%     drawnow
%     hold on
%     pause
end
%mse=sum((f-s).^2)/length(f);
prd=norm(f-s);%./sum((f-mean(f)).^2))*100;
mse=prd;

if show
    plot(s,'r','LineWidth',2);
    hold on;
    plot(f,'b','LineWidth',4);
    plot(s,'r','LineWidth',2);
    hold off;
end