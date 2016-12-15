%Legkisebb n�gyzetek m�dszere.%

%{
Param�terek:
 - x  : mintav�teli pontok x koordin�t�ja
 - y  : mintav�teli pontok y koordin�t�ja
 - order : illesztend� polinom foksz�ma
 - xx : hol akarjuk kisz�m�tani a polinom �rt�keit
Visszat�r�si �rt�kek:
 - s : az illesztett polinom �rt�kei az xx helyeken
 - p : az illesztett polinom egy�tthat�i
 
P�lda:    
    xx=-pi:0.1:pi;
    yy=cos(xx).*xx;
    leastsquare(xx,yy+rand(1,length(yy)),5,xx);
%}


function [s,p]=leastsquare(x,y,order,xx)
N=length(y);
if N<=order
    error('Helytelen param�terez�s(length(y)<order)!')
end
%�ltal�nos�tott inverz elk�sz�t�se%
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