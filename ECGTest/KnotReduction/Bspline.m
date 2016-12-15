%De Boor algoritmus l-ed fokú Bsplinok elõállítására%

%{
Paraméterek:
 - l  : fokszám
 - k  : megadja, hogy melyik alappontra illeszük a splinet
 - xk : alappontok x koordinátája
 - x  : ezeken a helyeken számoljuk ki a Bspline értékét
Visszatérési értékek:
 - y : a Bspline értéke az x helyeken
%}

function y=Bspline(l,k,xk,x)
if 0==l 
%    y=(x-xk(k))./(xk(k+1)-xk(k)) .* (x>=xk(k)).*(x<=xk(k+1));
%    y=y+(xk(k+2)-x)./(xk(k+2)-xk(k+1)) .* (x>xk(k+1)).*(x<=xk(k+2));
    y=ones(size(x)).*(x>=xk(k)).*(x<xk(k+1));
    
else
    %de Boor algoritmus%
    if (xk(k+l)-xk(k))
        y1=Bspline(l-1,k,xk,x).*(x-xk(k))./(xk(k+l)-xk(k));
    else
        y1=zeros(size(x));
    end
    if (xk(k+l+1)-xk(k+1))
        y2=Bspline(l-1,k+1,xk,x).*(xk(k+l+1)-x)./(xk(k+l+1)-xk(k+1));
    else
        y2=zeros(size(x));
    end
    y=y1 + y2;
end
