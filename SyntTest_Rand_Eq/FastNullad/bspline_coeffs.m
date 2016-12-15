%A n�gyzetesen legjobban k�zel�t� B-spline approxim�ci� meghat�roz�sa.

%{
Param�terek:
 - sig   : k�zel�tend� jel
 - knots : B-spline-ok alappontjai
 - order : B-spline-ok rendje (foksz�m + 1)
 - show  : logiai param�ter, ha igaz, akkor megjelen�ti az eredm�nyt
Visszat�r�si �rt�kek:
 - c     : n�gyzetesen legjobban k�zel�t� B-spline approxim�ci� egy�tthat�i
 - aprx  : n�gyzetesen legjobban k�zel�t� B-spline approxim�ci�
%}
function [c,aprx]=bspline_coeffs(sig,knots,order,show)
    x=1:1:length(sig);
    y=reshape(sig,length(sig),1);
    E=calcmat(order,x,knots);
    [Qll,Rll]=qr(E);
    zll=Qll'*y;
    c=Rll\zll;

    aprx=zeros(1,length(sig));
    %Peremfelt�tel: az els� �s utols� order-1 darab alappont egyezzen mem
    %az els� �s utols� alappontokkal. Az alappontok �j t�mbje 'tt' lesz.
    tt=[knots(1).*ones(order-1,1); knots; knots(end).*ones(order-1,1)]; 
    %Az els� �s utols� alappontban az approxim�ci�nak interpol�lnia kell.
    %Ez�rt a jelet �gy �ll�tjuk be/normaliz�ljuk, hogy az els� �s utols�
    %alappontban a f�ggv�ny �rt�ke 0 legyen. Ekkor az interpol�ci�s felt�telt 
    %figyelembe v�ve az els� �s utols� egy�tthat� 0. 
    c=[0; c; 0]; 
    for i=1:1:length(c)
        p=c(i)*BsplineDeBoor(order-1,i,tt,x);
        aprx=aprx+p;
    end
    
    %Kirajzol�s
    if show
        plot(aprx,'r','LineWidth',2);
        hold on;
        plot(sig,'b','LineWidth',4);
        plot(aprx,'r','LineWidth',2);
        hold off;
    end
end