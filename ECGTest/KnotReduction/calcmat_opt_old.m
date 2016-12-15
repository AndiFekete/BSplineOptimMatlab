%Optimalized version.%
function E=calcmat_opt_old(El,k,x,t,pp)
m=length(t)-2;
n=length(x)-2;
p=pp-1;
if p<1 || p>m
    error('The first and the last knot must be kept!');
end
%Discrete grid + boundary conditions%
tt=[t(1).*ones(1,k-1) t t(end).*ones(1,k-1)];

I1=ones(1,p-2);
I2=ones(1,m-p-1);
lambda=(t(pp)-tt(p:p+k-1)) ./ (tt(p+k:p+2*k-1)-tt(p:p+k-1));
mu=(tt(p+k+1:p+2*k)-t(pp)) ./ (tt(p+k+1:p+2*k)-tt(p+1:p+k));

i=zeros(1,m+2*k-3);
j=zeros(1,m+2*k-3);
s=zeros(1,m+2*k-3);
%Constructing the upper left corner (I1) of the 'B' matrix.%
i(1:p-2)=1:p-2;
j(1:p-2)=1:p-2;
s(1:p-2)=I1;

%B = sparse(1:p-2,1:p-2,I1,m+k-1,m+k-2);

%Constructing the middle part (B22) of the 'B' matrix.%
%B = sparse(p-1:p+k-2,p-1:p+k-2,lambda,m+k-1,m+k-2);
%B = sparse(p:p+k-1,p-1:p+k-2,mu,m+k-1,m+k-2);
i(p-1:p+k-2)=p-1:p+k-2;
j(p-1:p+k-2)=p-1:p+k-2;
s(p-1:p+k-2)=lambda;
% i(p:p+k-1)=p:p+k-1;
% j(p:p+k-1)=p:p+k-1;
% s(p:p+k-1)=lambda;


i(p+k-1:p+2*k-2)=p:p+k-1
j(p+k-1:p+2*k-2)=p-1:p+k-2
s(p+k-1:p+2*k-2)=mu;
% i(p+k:p+2*k-1)=p+1:p+k
% j(p+k:p+2*k-1)=p:p+k-1
% s(p+k:p+2*k-1)=mu;


%Constructing the bottom left corner (I2) of the 'B' matrix.%
%B = sparse(k+p+1:m+k-1,k+p:m+k-2,I2,m+k-1,m+k-2);
i(p+2*k-1:end)=k+p:m+k-2;
j(p+2*k-1:end)=k+p-1:m+k-3;
s(p+2*k-1:end)=I2;

B = sparse(i,j,s,m+k-2,m+k-3)

spy(B)
E=El*B;

% E=zeros(n,m+k-2);
% for i=1:(m+k-2)
%     E(:,i)=Bspline(k-1,i+1,tt,x(2:end-1))';
% %    plot(linspace(x(1),x(end),100),Bspline(k-1,i+1,tt,linspace(x(1),x(end),100)))
% %     drawnow
% %     hold on
% %     pause
% end