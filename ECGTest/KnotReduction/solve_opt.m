function [c,Rll,zll]=solve_opt(Rl,zl,B,b,k,p)
Rll=Rl*B;
zll=zl;
[m,n]=size(Rll);


if p-2<=0
    ind=1;
else
    ind=p-2;
end

for i=ind:1:n-k-1
    [G,y] = planerot(Rll(i:i+1,i));
    %Rotating vector 'zll'.%
    r11=G(1,1)*zll(i)+G(1,2)*zll(i+1);
    r12=G(2,1)*zll(i)+G(2,2)*zll(i+1);    
    zll(i:i+1)=[r11,r12];    
    %Rotating matrix 'Rll'.%
    Rll(i:i+1,i)=y;        
    for j=1:1:k-1
        r11=G(1,1)*Rll(i,i+j)+G(1,2)*Rll(i+1,i+j);
        r12=G(2,1)*Rll(i,i+j)+G(2,2)*Rll(i+1,i+j);    
        Rll(i:i+1,i+j)=[r11,r12];
    end
end    

%if i>n-k%
if n-k<=0
    ind=1;
else
    ind=n-k;
end

for i=ind:1:n
    [G,y] = planerot(Rll(i:i+1,i));
    %Rotating vector 'zll'.%
    zll(i:i+1)=G*zll(i:i+1);    
    %Rotating matrix 'Rll'.%
    Rll(i:i+1,i)=y;        
    for j=1:1:n-i
        Rll(i:i+1,i+j)=G*Rll(i:i+1,i+j);
    end
end    

%c=Rll\zll;
%c=Rll(1:n,1:n)\zll(1:n);
opt.UT=true;
c = linsolve(Rll(1:n,1:n),zll(1:n),opt);

