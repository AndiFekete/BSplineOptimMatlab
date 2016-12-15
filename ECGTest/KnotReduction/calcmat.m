function E=calcmat(k,x,t)
%Discrete grid + boundary conditions%
tt=[t(1).*ones(1,k-1) t t(end).*ones(1,k-1)];
n=length(x)-2;
m=length(t)-2;

E=zeros(n,m+k-2);
for i=1:(m+k-2)
    E(:,i)=Bspline(k-1,i+1,tt,x(2:end-1))';
%    plot(linspace(x(1),x(end),100),Bspline(k-1,i+1,tt,linspace(x(1),x(end),100)))
%     drawnow
%     hold on
%     pause
end