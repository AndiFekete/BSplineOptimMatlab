%K�b�s splinok el��ll�t�sa Bspline b�zisban%

%{
Param�terek:
 - xk : alappontok x koordin�t�ja
 - yk : alappontok y koordin�t�ja
 - Da : els� v�gpontban a deriv�lt �rt�ke
 - Db : m�sodik v�gpontban a deriv�lt �rt�ke
 - x  : ezeken a helyeken szeretn�nk kisz�m�tani a spline �rt�k�t
Visszat�r�si �rt�kek:
 - y : a Bspline b�zisban fel�rt spline �rt�ke az x helyeken

P�lda:
  f=@(x) cos(x).*x;
  Df=@(x) -sin(x).*x+cos(x);
  x=linspace(1,5,100);
  xk=[ 1 2 3 4 5]
  yk=cos(xk).*xk;
  Bspline3(xk,yk,Df(xk(1)),f(xk(end)),x)
%}
function y=Bspline3(xk,yk,Da,Db,x)
h=xk(2)-xk(1);
n=length(xk);
d=[0 4.*ones(1,n) 0];
l=[ones(1,n) 0];
u=[0 ones(1,n)];
i=[1:n+2 2:n+2 1:n+1];
j=[1:n+2 1:n+1 2:n+2];
A=sparse(i,j,[d,l,u]);
%Hermite-f�le peremfelt�tel hozz�ad�sa az A m�trixhoz%
A(1,1:3)=[-3/h 0 3/h];
A(end,end-2:end)=[-3/h 0 3/h];
b=6*[Da yk Db];
%spy(A);
full(A);
c=A\b';
S=zeros(1,length(x));
Xk=[xk(1)-3*h xk(1)-2*h xk(1)-h xk xk(end)+h xk(end)+2*h xk(end)+3*h];
for i=1:1:n+2
    S=S+c(i).*Bspline(3,i,Xk,x);
end
plot(x,S)
hold on m
plot(xk,yk,'r.')