%l-ed fok� Bsplinok el��ll�t�sa osztott differenci�k seg�ts�g�vel.%

%{
Param�terek:
 - l  : foksz�m
 - k  : megadja, hogy melyik alappontra illesz�k a splinet
 - xk : alappontok x koordin�t�ja
 - x  : ezeken a helyeken sz�moljuk ki a Bspline �rt�k�t
Visszat�r�si �rt�kek:
 - y : a Bspline �rt�ke az x helyeken
%}

function y=BsplineDeBoor(l,k,xk,x)
    tpvals=trunpow(l,xk(k:k+l+1),x);
    y=tpvals;
    for i=1:1:l+1
        y=divdiff(l,i,y,x,xk(k:k+l+1));
    end
    y=(-1)^(l+1)*y; %scaling
    %y(x<xk(1) | xk(end)<x)=0;
end

%Trucated power functions.
function val=trunpow(l,xj,x)
    val=zeros(length(xj),length(x));
    for j=1:1:length(xj)
        y=(x-xj(j));
        val(j,y>0)=y(y>0).^l;
    end
end

%Kth derivative of the trucated power functions.
function val=dertrunpow(l,k,xj,x)
    val=zeros(length(xj),length(x));
    for j=1:1:length(xj)
        y=(x-xj(j));
        val(j,y>0)=prod((l-k+1):l)*y(y>0).^(l-k)*(-1)^k;
    end
end

%Kth order divided differences between rows.
function d=divdiff(l,k,y,x,xj)
    d=diff(y,1,1);
    for i=1:1:size(d,1)
        if xj(i+k)==xj(i)
            d(i,:)=dertrunpow(l,k,xj(i),x)./factorial(k);
        else
            h=xj(i+k)-xj(i);
            d(i,:)=d(i,:)/h;
        end
    end
end