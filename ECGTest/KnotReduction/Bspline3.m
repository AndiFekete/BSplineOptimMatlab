%Köbös splinok elõállítása Bspline bázisban%

%{
Paraméterek:
 - xk : alappontok x koordinátája
 - yk : alappontok y koordinátája
 - Da : elsõ végpontban a derivált értéke
 - Db : második végpontban a derivált értéke
 - x  : ezeken a helyeken szeretnénk kiszámítani a spline értékét
Visszatérési értékek:
 - y : a Bspline bázisban felírt spline értéke az x helyeken

Példa:
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
%Hermite-féle peremfeltétel hozzáadása az A mátrixhoz%
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