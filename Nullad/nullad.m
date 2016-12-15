function [knots,aprx,err]=nullad(x,y,k,show)
    
    %Initialize variables
    N=length(x);
    F=cumsum(y);
    l=1;  %initial number of subintervals
     
    %Initialize the boundary knots
    knots=[x(2), x(N-1)];
    ind=[2,N-1];
    
    %% Effective implementation of the knot ranking algorithm for zero order B-spline approximation
    for i=1:1:k
        maxok=zeros(1,l);
        maxind=zeros(1,l);
        Fdiff=F(ind(2:end)-1)-F(ind(1:end-1)-1);
        wavgFdiff=sum(Fdiff.^2 ./ (knots(2:end)-knots(1:end-1)));  %weighted average of the square of the subintegrals 
        erreffective=zeros(1,N);
        for j=1:l
                actFdiff=(F(ind(j+1)-1)-F(ind(j)-1));              %Fdiff of the actually "divided" subinterval
                actwavgFdiff=actFdiff.^2 ./ (knots(j+1)-knots(j)); %weighted average of the square of the actual subintegral 
                aprxerrj=wavgFdiff-actwavgFdiff;                   %approximation error outside the actually divided subinterval
                p1=ind(j);
                p2=ind(j+1);
                
                if p2-p1>1  %if the subinterval contains more than 1 samples
                    values=r(F,x,p1+1:p2-1,p1,p2,actFdiff,l); %ha ind(j),ind(j+1) közé nem szúrtunk be alappontot 
                                                              %akkor ez ugyan az mint az elõzõ iterációban, ki sem kellene számolni           
                    erreffective(p1+1:p2-1)=values+aprxerrj;
                    [maxok(j),maxind(j)]=max(values+aprxerrj);
                    maxind(j)=maxind(j)+p1;     
                end
        end
        
        [m,mi]=max(maxok);
        maxindex = maxind(mi);
        elem=x(maxindex);
        knots=sort([knots elem]);
        ind=sort([ind maxindex]);
        l=l+1;
        %% Displaying the approximation
        if show
            [c,err,aprx]=bspline_coeffs(y,[x(1) knots x(N)] ,1);
            plot(x,y,'b',x,aprx,'r');
            legend('Original',sprintf('Approx. (PRD: %.1f%%)',err));
            drawnow;
        end        
    end
    [c,err,aprx]=bspline_coeffs(y,[x(1) knots x(N)] ,1);    
end

%% Auxiliary functions

function v=r(F,x,i,p1,p2,norma,l)
    r1=(F(i-1)-F(p1-1)).^2./(x(i)-x(p1));
    r2=(norma-(F(i-1)-F(p1-1))).^2./(x(p2)-x(i));    
    v=r1+r2;
end
    
    
 