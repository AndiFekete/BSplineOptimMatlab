%2-od fok� Bsplinok alappont szerinti deriv�ltj�nak el��ll�t�sa explicit alakban%

%{
Param�terek:
 - k  : megadja, hogy melyik alappontra illesz�k a splinet
 - j  : megadja, hogy hanyadik alappont szerint deriv�ljunk
 - xk : alappontok x koordin�t�ja
 - x  : ezeken a helyeken sz�moljuk ki a Bspline �rt�k�t
Visszat�r�si �rt�kek:
 - y : a Bspline �rt�ke az x helyeken
%}

function y=DerBspline2(k,j,xk,x)
    switch j-k 
        case 0 % xk(k) alappont szerinti deriv�lt
            xx=x(xk(k)<=x & x<xk(k+1));     %A felvett �rt�kek az [xk(k),xk(k+1)] r�szintervallumon.
            d1=(xk(k+1).*xk(k+2)-xk(k).*xk(k+2)-xk(k).*xk(k+1) + xk(k).^2);
            y1=(2*(xk(k)-xx)*d1 - (xx.^2 - 2*xx.*xk(k) + xk(k).^2).*(2*xk(k)-xk(k+2)-xk(k+1))) / d1^2;

            xx=x(xk(k+1)<=x & x<xk(k+2));   %A felvett �rt�kek az [xk(k+1),xk(k+2)] r�szintervallumon.
            d2=(xk(k).*xk(k+1)-xk(k).*xk(k+2)-xk(k+1).*xk(k+2) + xk(k+2).^2);
            y2=((xx-xk(k+2))*d2 - (xk(k+2).*xx + xk(k).*xx - xk(k).*xk(k+2) - xx.^2).*(xk(k+1)-xk(k+2))) / d2^2;

            y=zeros(size(x));
            y(xk(k)<=x & x<xk(k+2))=[y1,y2];
        case 1 % xk(k+1) alappont szerinti deriv�lt
            xx=x(xk(k)<=x & x<xk(k+1));     %A felvett �rt�kek az [xk(k),xk(k+1)] r�szintervallumon.
            d1=(xk(k+1).*xk(k+2)-xk(k).*xk(k+2)-xk(k).*xk(k+1) + xk(k).^2);
            y1=(-(xx.^2 - 2*xx.*xk(k) + xk(k).^2).*(xk(k+2) - xk(k))) / d1^2;

            xx=x(xk(k+1)<=x & x<xk(k+2));   %A felvett �rt�kek az [xk(k+1),xk(k+2)] r�szintervallumon.
            d2=(xk(k).*xk(k+1)-xk(k).*xk(k+2)-xk(k+1).*xk(k+2) + xk(k+2).^2);
            d3=(xk(k+2).*xk(k+3)-xk(k+1).*xk(k+3)-xk(k+1).*xk(k+2) + xk(k+1).^2);
            y2=-(xk(k+2).*xx + xk(k).*xx - xk(k).*xk(k+2) - xx.^2).*(xk(k) - xk(k+2)) / d2^2 + ...
               +((xx - xk(k+3)).*d3) / d3^2 - (xk(k+3).*xx + xk(k+1).*xx - xk(k+1).*xk(k+3) - xx.^2).*(2*xk(k+1) - xk(k+3) - xk(k+2)) / d3^2;  

            xx=x(xk(k+2)<=x & x<xk(k+3));   %A felvett �rt�kek az [xk(k+2),xk(k+3)] r�szintervallumon.
            d4=(xk(k+1).*xk(k+2)-xk(k+1).*xk(k+3)-xk(k+2).*xk(k+3) + xk(k+3).^2);
            y3=(-(xk(k+3).^2 - 2.*xx.*xk(k+3) + xx.^2).*(xk(k+2) - xk(k+3))) / d4^2;
           
            y=zeros(size(x));
            y(xk(k)<=x & x<xk(k+3))=[y1,y2,y3];

        case 2 % xk(k+2) alappont szerinti deriv�lt
            xx=x(xk(k)<=x & x<xk(k+1));     %A felvett �rt�kek az [xk(k),xk(k+1)] r�szintervallumon.
            d1=(xk(k+1).*xk(k+2)-xk(k).*xk(k+2)-xk(k).*xk(k+1) + xk(k).^2);
            y1=(-(xx.^2 - 2*xx.*xk(k) + xk(k).^2).*(xk(k+1) - xk(k))) / d1^2;

            xx=x(xk(k+1)<=x & x<xk(k+2));   %A felvett �rt�kek az [xk(k+1),xk(k+2)] r�szintervallumon.
            d2=(xk(k).*xk(k+1)-xk(k).*xk(k+2)-xk(k+1).*xk(k+2) + xk(k+2).^2);
            d3=(xk(k+2).*xk(k+3)-xk(k+1).*xk(k+3)-xk(k+1).*xk(k+2) + xk(k+1).^2);
            y2=((xx-xk(k))*d2 - (xk(k+2).*xx + xk(k).*xx - xk(k)*xk(k+2) - xx.^2).*(2*xk(k+2) - xk(k) - xk(k+1))) ./ d2^2 + ...
                - (xk(k+3).*xx + xk(k+1).*xx - xk(k+1).*xk(k+3) - xx.^2).*(xk(k+3)-xk(k+1)) ./ d3^2;  

            xx=x(xk(k+2)<=x & x<xk(k+3));   %A felvett �rt�kek az [xk(k+2),xk(k+3)] r�szintervallumon.
            d4=(xk(k+1).*xk(k+2)-xk(k+1).*xk(k+3)-xk(k+2).*xk(k+3) + xk(k+3).^2);
            y3=(-(xk(k+3).^2 - 2.*xx.*xk(k+3) + xx.^2).*(xk(k+1) - xk(k+3))) / d4^2;
           
            y=zeros(size(x));
            y(xk(k)<=x & x<xk(k+3))=[y1,y2,y3];
            
            
        case 3 % xk(k+3) alappont szerinti deriv�lt
            xx=x(xk(k+1)<=x & x<xk(k+2));     %A felvett �rt�kek az [xk(k+1),xk(k+2)] r�szintervallumon.
            d1=(xk(k+2).*xk(k+3)-xk(k+1).*xk(k+3)-xk(k+1).*xk(k+2) + xk(k+1).^2);
            y1=((xx-xk(k+1))*d1 - (xk(k+3).*xx + xk(k+1).*xx - xk(k+1).*xk(k+3) - xx.^2).*(xk(k+2)-xk(k+1))) / d1^2;

            xx=x(xk(k+2)<=x & x<xk(k+3));   %A felvett �rt�kek az [xk(k+2),xk(k+3)] r�szintervallumon.
            d2=(xk(k+1).*xk(k+2)-xk(k+1).*xk(k+3)-xk(k+2).*xk(k+3) + xk(k+3).^2);
            y2=(2*(xk(k+3)-xx)*d2 - (xk(k+3).^2 - 2*xx.*xk(k+3) + xx.^2).*(2*xk(k+3)-xk(k+1)-xk(k+2))) / d2^2;

            y=zeros(size(x));
            y(xk(k+1)<=x & x<xk(k+3))=[y1,y2];            
        otherwise
            y=zeros(size(x));
    end        
end
