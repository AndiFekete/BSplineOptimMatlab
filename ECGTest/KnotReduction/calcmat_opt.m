%Optimalized version.%
function [B,E]=calcmat_opt2(k,x,t,pp,El)
m=length(t)-2;
n=length(x)-2;
p=pp-1;
if p<1 || p>m
    error('The first and the last knot must be kept!');
end

%Discrete grid + boundary conditions%
tt=[t(1).*ones(1,k-1) t t(end).*ones(1,k-1)];

I1=ones(1,p-1);
I2=ones(1,m-p);
lambda=(t(pp)-tt(p:p+k-1)) ./ (tt(p+k:p+2*k-1)-tt(p:p+k-1));
mu=(tt(p+k+1:p+2*k)-t(pp)) ./ (tt(p+k+1:p+2*k)-tt(p+1:p+k));
i=zeros(1,m+2*k-1);
j=zeros(1,m+2*k-1);
s=zeros(1,m+2*k-1);

%Constructing the upper left corner (I1) of the 'B' matrix.%
i(1:p-1)=1:p-1;
j(1:p-1)=1:p-1;
s(1:p-1)=I1;

%Constructing the middle part (B22) of the 'B' matrix.%
i(p:p+k-1)=p:p+k-1;
j(p:p+k-1)=p:p+k-1;
s(p:p+k-1)=lambda;

i(p+k:p+2*k-1)=p+1:p+k;
j(p+k:p+2*k-1)=p:p+k-1;
s(p+k:p+2*k-1)=mu;


%Constructing the bottom left corner (I2) of the 'B' matrix.%
i(p+2*k:end)=k+p+1:m+k;
j(p+2*k:end)=k+p:m+k-1;
s(p+2*k:end)=I2;

%Deleting unnecessary rows and columns.%
B = sparse(i,j,s,m+k,m+k-1);
B(1,:)=[];
B(:,1)=[];
B(end,:)=[];
B(:,end)=[];

if nargin>4
    E=El*B;
end