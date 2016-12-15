%2-od fokú Bsplinok elõállítása explicit alakban%

%{
Paraméterek:
 - k  : megadja, hogy melyik alappontra illeszük a splinet
 - xk : alappontok x koordinátája
 - x  : ezeken a helyeken számoljuk ki a Bspline értékét
Visszatérési értékek:
 - y : a Bspline értéke az x helyeken
%}

function y=Bspline2(k,xk,x)
    xx=x(xk(k)<=x & x<xk(k+1));   %A felvett értékek az [xk(k),xk(k+1)] részintervallumon.
    y1=(xx.^2 - 2*xx.*xk(k) + xk(k).^2) ./ (xk(k+1).*xk(k+2)-xk(k).*xk(k+2)-xk(k).*xk(k+1) + xk(k).^2); 

    xx=x(xk(k+1)<=x & x<xk(k+2)); %A felvett értékek az [xk(k+1),xk(k+2)] részintervallumon.
    a=(xk(k+2).*xx + xk(k).*xx - xk(k).*xk(k+2) - xx.^2) ./ (xk(k).*xk(k+1)-xk(k).*xk(k+2)-xk(k+1).*xk(k+2) + xk(k+2).^2); 
    b=(xk(k+3).*xx + xk(k+1).*xx - xk(k+1).*xk(k+3) - xx.^2) ./ (xk(k+2).*xk(k+3)-xk(k+1).*xk(k+3)-xk(k+1).*xk(k+2) + xk(k+1).^2); 
    y2=a+b;
    
    xx=x(xk(k+2)<=x & x<xk(k+3)); %A felvett értékek az [xk(k+2),xk(k+3)] részintervallumon.
    y3=(xk(k+3).^2 - 2*xx.*xk(k+3) + xx.^2) ./ (xk(k+1).*xk(k+2)-xk(k+1).*xk(k+3)-xk(k+2).*xk(k+3) + xk(k+3).^2); 

    y=zeros(size(x));
    y(xk(k)<=x & x<xk(k+3))=[y1,y2,y3];
end
