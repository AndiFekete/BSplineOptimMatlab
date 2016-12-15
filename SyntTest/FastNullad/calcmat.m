%A t�lhat�rozott egyenletrendszer m�trixa n�gyzetesen legjobban k�zel�t� B-spline-okhoz.

%{
Param�terek:
 - k  : B-spline-ok rendje (foksz�m + 1)
 - x  : ezeken a helyeken sz�moljuk ki a Bspline �rt�k�t
 - t  : a B-spline-ok alappontjai
Visszat�r�si �rt�kek:
 - E : az egyenletrendszer m�trixa (E_ij= j. B-spline �rt�ke a i. mint�nak megfelel� ponton)
%}
function E=calcmat(k,x,t)
%Discrete grid + boundary conditions%
tt=[t(1).*ones(k-1,1); t; t(end).*ones(k-1,1)];
n=length(x);
m=length(t)-2;

E=zeros(n,m+k-2);
for i=1:(m+k-2)
    E(:,i)=BsplineDeBoor(k-1,i+1,tt,x)';
%     plot(linspace(x(1),x(end),100),Bspline(k-1,i+1,tt,linspace(x(1),x(end),100)),'b',tt,Bspline(k-1,i+1,tt,tt),'r.');    
%     drawnow
%     hold on;
%     pause
end