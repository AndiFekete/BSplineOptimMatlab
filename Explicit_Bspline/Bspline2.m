%2-od fok� Bsplinok el��ll�t�sa explicit alakban%

%{
Param�terek:
 - k  : megadja, hogy melyik alappontra illesz�k a splinet
 - xk : alappontok x koordin�t�ja
 - x  : ezeken a helyeken sz�moljuk ki a Bspline �rt�k�t
Visszat�r�si �rt�kek:
 - y : a Bspline �rt�ke az x helyeken
%}

function y=Bspline2(k,xk,x)
    xx=x(xk(k)<=x & x<xk(k+1));   %A felvett �rt�kek az [xk(k),xk(k+1)] r�szintervallumon.
    y1=(xx.^2 - 2*xx.*xk(k) + xk(k).^2) ./ (xk(k+1).*xk(k+2)-xk(k).*xk(k+2)-xk(k).*xk(k+1) + xk(k).^2); 

    xx=x(xk(k+1)<=x & x<xk(k+2)); %A felvett �rt�kek az [xk(k+1),xk(k+2)] r�szintervallumon.
    a=(xk(k+2).*xx + xk(k).*xx - xk(k).*xk(k+2) - xx.^2) ./ (xk(k).*xk(k+1)-xk(k).*xk(k+2)-xk(k+1).*xk(k+2) + xk(k+2).^2); 
    b=(xk(k+3).*xx + xk(k+1).*xx - xk(k+1).*xk(k+3) - xx.^2) ./ (xk(k+2).*xk(k+3)-xk(k+1).*xk(k+3)-xk(k+1).*xk(k+2) + xk(k+1).^2); 
    y2=a+b;
    
    xx=x(xk(k+2)<=x & x<xk(k+3)); %A felvett �rt�kek az [xk(k+2),xk(k+3)] r�szintervallumon.
    y3=(xk(k+3).^2 - 2*xx.*xk(k+3) + xx.^2) ./ (xk(k+1).*xk(k+2)-xk(k+1).*xk(k+3)-xk(k+2).*xk(k+3) + xk(k+3).^2); 

    y=zeros(size(x));
    y(xk(k)<=x & x<xk(k+3))=[y1,y2,y3];
end
