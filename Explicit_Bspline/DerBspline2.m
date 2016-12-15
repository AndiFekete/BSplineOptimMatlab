%2-od fokú Bsplinok alappont szerinti deriváltjának elõállítása explicit alakban%

%{
Paraméterek:
 - k  : megadja, hogy melyik alappontra illeszük a splinet
 - j  : megadja, hogy hanyadik alappont szerint deriváljunk
 - xk : alappontok x koordinátája
 - x  : ezeken a helyeken számoljuk ki a Bspline értékét
Visszatérési értékek:
 - y : a Bspline értéke az x helyeken
%}

function y=DerBspline2(k,j,xk,x)
    switch j-k 
        case 0 % xk(k) alappont szerinti derivált
            xx=x(xk(k)<=x & x<xk(k+1));     %A felvett értékek az [xk(k),xk(k+1)] részintervallumon.
            d1=(xk(k+1).*xk(k+2)-xk(k).*xk(k+2)-xk(k).*xk(k+1) + xk(k).^2);
            y1=(2*(xk(k)-xx)*d1 - (xx.^2 - 2*xx.*xk(k) + xk(k).^2).*(2*xk(k)-xk(k+2)-xk(k+1))) / d1^2;

            xx=x(xk(k+1)<=x & x<xk(k+2));   %A felvett értékek az [xk(k+1),xk(k+2)] részintervallumon.
            d2=(xk(k).*xk(k+1)-xk(k).*xk(k+2)-xk(k+1).*xk(k+2) + xk(k+2).^2);
            y2=((xx-xk(k+2))*d2 - (xk(k+2).*xx + xk(k).*xx - xk(k).*xk(k+2) - xx.^2).*(xk(k+1)-xk(k+2))) / d2^2;

            y=zeros(size(x));
            y(xk(k)<=x & x<xk(k+2))=[y1,y2];
        case 1 % xk(k+1) alappont szerinti derivált
            xx=x(xk(k)<=x & x<xk(k+1));     %A felvett értékek az [xk(k),xk(k+1)] részintervallumon.
            d1=(xk(k+1).*xk(k+2)-xk(k).*xk(k+2)-xk(k).*xk(k+1) + xk(k).^2);
            y1=(-(xx.^2 - 2*xx.*xk(k) + xk(k).^2).*(xk(k+2) - xk(k))) / d1^2;

            xx=x(xk(k+1)<=x & x<xk(k+2));   %A felvett értékek az [xk(k+1),xk(k+2)] részintervallumon.
            d2=(xk(k).*xk(k+1)-xk(k).*xk(k+2)-xk(k+1).*xk(k+2) + xk(k+2).^2);
            d3=(xk(k+2).*xk(k+3)-xk(k+1).*xk(k+3)-xk(k+1).*xk(k+2) + xk(k+1).^2);
            y2=-(xk(k+2).*xx + xk(k).*xx - xk(k).*xk(k+2) - xx.^2).*(xk(k) - xk(k+2)) / d2^2 + ...
               +((xx - xk(k+3)).*d3) / d3^2 - (xk(k+3).*xx + xk(k+1).*xx - xk(k+1).*xk(k+3) - xx.^2).*(2*xk(k+1) - xk(k+3) - xk(k+2)) / d3^2;  

            xx=x(xk(k+2)<=x & x<xk(k+3));   %A felvett értékek az [xk(k+2),xk(k+3)] részintervallumon.
            d4=(xk(k+1).*xk(k+2)-xk(k+1).*xk(k+3)-xk(k+2).*xk(k+3) + xk(k+3).^2);
            y3=(-(xk(k+3).^2 - 2.*xx.*xk(k+3) + xx.^2).*(xk(k+2) - xk(k+3))) / d4^2;
           
            y=zeros(size(x));
            y(xk(k)<=x & x<xk(k+3))=[y1,y2,y3];

        case 2 % xk(k+2) alappont szerinti derivált
            xx=x(xk(k)<=x & x<xk(k+1));     %A felvett értékek az [xk(k),xk(k+1)] részintervallumon.
            d1=(xk(k+1).*xk(k+2)-xk(k).*xk(k+2)-xk(k).*xk(k+1) + xk(k).^2);
            y1=(-(xx.^2 - 2*xx.*xk(k) + xk(k).^2).*(xk(k+1) - xk(k))) / d1^2;

            xx=x(xk(k+1)<=x & x<xk(k+2));   %A felvett értékek az [xk(k+1),xk(k+2)] részintervallumon.
            d2=(xk(k).*xk(k+1)-xk(k).*xk(k+2)-xk(k+1).*xk(k+2) + xk(k+2).^2);
            d3=(xk(k+2).*xk(k+3)-xk(k+1).*xk(k+3)-xk(k+1).*xk(k+2) + xk(k+1).^2);
            y2=((xx-xk(k))*d2 - (xk(k+2).*xx + xk(k).*xx - xk(k)*xk(k+2) - xx.^2).*(2*xk(k+2) - xk(k) - xk(k+1))) ./ d2^2 + ...
                - (xk(k+3).*xx + xk(k+1).*xx - xk(k+1).*xk(k+3) - xx.^2).*(xk(k+3)-xk(k+1)) ./ d3^2;  

            xx=x(xk(k+2)<=x & x<xk(k+3));   %A felvett értékek az [xk(k+2),xk(k+3)] részintervallumon.
            d4=(xk(k+1).*xk(k+2)-xk(k+1).*xk(k+3)-xk(k+2).*xk(k+3) + xk(k+3).^2);
            y3=(-(xk(k+3).^2 - 2.*xx.*xk(k+3) + xx.^2).*(xk(k+1) - xk(k+3))) / d4^2;
           
            y=zeros(size(x));
            y(xk(k)<=x & x<xk(k+3))=[y1,y2,y3];
            
            
        case 3 % xk(k+3) alappont szerinti derivált
            xx=x(xk(k+1)<=x & x<xk(k+2));     %A felvett értékek az [xk(k+1),xk(k+2)] részintervallumon.
            d1=(xk(k+2).*xk(k+3)-xk(k+1).*xk(k+3)-xk(k+1).*xk(k+2) + xk(k+1).^2);
            y1=((xx-xk(k+1))*d1 - (xk(k+3).*xx + xk(k+1).*xx - xk(k+1).*xk(k+3) - xx.^2).*(xk(k+2)-xk(k+1))) / d1^2;

            xx=x(xk(k+2)<=x & x<xk(k+3));   %A felvett értékek az [xk(k+2),xk(k+3)] részintervallumon.
            d2=(xk(k+1).*xk(k+2)-xk(k+1).*xk(k+3)-xk(k+2).*xk(k+3) + xk(k+3).^2);
            y2=(2*(xk(k+3)-xx)*d2 - (xk(k+3).^2 - 2*xx.*xk(k+3) + xx.^2).*(2*xk(k+3)-xk(k+1)-xk(k+2))) / d2^2;

            y=zeros(size(x));
            y(xk(k+1)<=x & x<xk(k+3))=[y1,y2];            
        otherwise
            y=zeros(size(x));
    end        
end
