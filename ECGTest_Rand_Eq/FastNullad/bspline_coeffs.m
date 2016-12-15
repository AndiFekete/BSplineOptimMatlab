%A négyzetesen legjobban közelítõ B-spline approximáció meghatározása.

%{
Paraméterek:
 - sig   : közelítendõ jel
 - knots : B-spline-ok alappontjai
 - order : B-spline-ok rendje (fokszám + 1)
 - show  : logiai paraméter, ha igaz, akkor megjeleníti az eredményt
Visszatérési értékek:
 - c     : négyzetesen legjobban közelítõ B-spline approximáció együtthatói
 - aprx  : négyzetesen legjobban közelítõ B-spline approximáció
%}
function [c,aprx]=bspline_coeffs(sig,knots,order,show)
    x=1:1:length(sig);
    y=reshape(sig,length(sig),1);
    E=calcmat(order,x,knots);
    [Qll,Rll]=qr(E);
    zll=Qll'*y;
    c=Rll\zll;

    aprx=zeros(1,length(sig));
    %Peremfeltétel: az elsõ és utolsó order-1 darab alappont egyezzen mem
    %az elsõ és utolsó alappontokkal. Az alappontok új tömbje 'tt' lesz.
    tt=[knots(1).*ones(order-1,1); knots; knots(end).*ones(order-1,1)]; 
    %Az elsõ és utolsó alappontban az approximációnak interpolálnia kell.
    %Ezért a jelet úgy állítjuk be/normalizáljuk, hogy az elsõ és utolsó
    %alappontban a függvény értéke 0 legyen. Ekkor az interpolációs feltételt 
    %figyelembe véve az elsõ és utolsó együttható 0. 
    c=[0; c; 0]; 
    for i=1:1:length(c)
        p=c(i)*BsplineDeBoor(order-1,i,tt,x);
        aprx=aprx+p;
    end
    
    %Kirajzolás
    if show
        plot(aprx,'r','LineWidth',2);
        hold on;
        plot(sig,'b','LineWidth',4);
        plot(aprx,'r','LineWidth',2);
        hold off;
    end
end