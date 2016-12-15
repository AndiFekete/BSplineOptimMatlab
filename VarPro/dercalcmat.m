%A túlhatározott egyenletrendszer mátrixa négyzetesen legjobban közelítõ B-spline-okhoz.

%{
Paraméterek:
 - k       : B-spline-ok rendje (fokszám + 1)
 - x       : ezeken a helyeken számoljuk ki a Bspline értékét
 - t       : a B-spline-ok alappontjai
 - sortind : elõfordulhat, hogy az alappontok sorrendje megváltozik egy Gauss-Newton lépés után.
             Azt, hogy a rendezés után hova kerülnek az egyes változók a 'sortind' indextömb rögzíti.
Visszatérési értékek:
 - dPhi : a Phi_j bázisfüggvények t_i szerinti parciális parciális 
          deriváltfüggvényei, ahol j=Ind(1,k) és i=Ind(2,k).
%}
function [dPhi, Ind]=dercalcmat(k,x,t,sortind)
    %Discrete grid + boundary conditions%
    tt=[t(1).*ones(k-1,1); t; t(end).*ones(k-1,1)];
    n=length(x);
    m=length(t)-2;
    p1=((m+k-4)-k-2)*(k+1); %Number of partial derivatives for Bsplines with interiror points. 
    p2=2*(k*(k+1)/2);     %Number of partial derivatives for Bsplines which has boundary points. 
    p=p1+p2;   %Number of partial derivatives with and without (left & right) boundary knots.
    dPhi=zeros(n,p);
    Ind=zeros(2,p);

    %Partial derivatives of the Bsplines for which the left boundary point is a knot. 
    top=1;
    for j=2:k
        for i=1:j
            Ind(:,top)=[j-1;sortind(i)];
            dPhi(:,top)=DerBsplineDeBoor(k-1,j,k+i,tt,x)';
            top=top+1;
        end
    end

    %Partial derivatives of the Bsplines without boundary point knot.
    for j=k+1:(m+k-4)    
        for i=0:k
            Ind(:,top)=[j-1;sortind(j-k+i)];
            dPhi(:,top)=DerBsplineDeBoor(k-1,j,j+i,tt,x)';
            top=top+1;
        end
    end

    %Partial derivatives of the last Bsplines for which the right boundary point is a knot. 
    for j=m+1:m+k-1
        for i=m+k-j:-1:0
            Ind(:,top)=[j-1;sortind(j+i-k)];
            dPhi(:,top)=DerBsplineDeBoor(k-1,j,j+i,tt,x)';
            top=top+1;
        end
    end
end