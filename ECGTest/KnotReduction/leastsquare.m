%Legkisebb négyzetek módszere.%

%{
Paraméterek:
 - x  : mintavételi pontok x koordinátája
 - y  : mintavételi pontok y koordinátája
 - order : illesztendõ polinom fokszáma
 - xx : hol akarjuk kiszámítani a polinom értékeit
Visszatérési értékek:
 - s : az illesztett polinom értékei az xx helyeken
 - p : az illesztett polinom együtthatói
 
Példa:    
    xx=-pi:0.1:pi;
    yy=cos(xx).*xx;
    leastsquare(xx,yy+rand(1,length(yy)),5,xx);
%}


function [s,p]=leastsquare(x,y,order,xx)
N=length(y);
if N<=order
    error('Helytelen paraméterezés(length(y)<order)!')
end
%Általánosított inverz elkészítése%
A=zeros(N,order+1);
for i=1:1:order+1
    A(:,i)=x.^(i-1);
end
A
b=y;
Ap=inv(A'*A)*A'
p=Ap*b'
s=polyval(p(end:-1:1),xx);
plot(xx,s);
hold on;
plot(x,y,'r.');