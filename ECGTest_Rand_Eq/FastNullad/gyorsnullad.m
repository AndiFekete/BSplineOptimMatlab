function [knots,aprx,err]=gyorsnullad(x,y,k,show)
    N=length(x);
    
    %Initialize the boundary knots
    knots=[x(2), x(N-1)];
    ind=[2, N-1];

    %Initialize variables
    F=cumsum(y);
    maxok=zeros(1,k);
    maxind=zeros(1,k);
    wavgFdiffj=zeros(1,k);
    erreffective=zeros(1,N);        
    
    chosen=1;   %index of the chosen subinterval

    firstknot = nullad(x,y,1,0);           %choosing the first knot
    knots=sort([knots x(firstknot(2))]);   %inserting the first knot
    ind=sort([ind firstknot(2)]);          %inserting the first knot index
    
    %% Effective implementation for more than 2 knots
    for i=1:1:k-1    
        Fdiff=F(ind(2:end)-1)-F(ind(1:end-1)-1);
        chosenmax=[0 0];
        chosenind=[0 0];     
        chosenwavgFdiff=[0 0];
        for j=chosen:chosen+1
                actFdiff=(F(ind(j+1)-1)-F(ind(j)-1)); %Fdiff of the actually "divided" subinterval
                chosenwavgFdiff(j-chosen+1)=actFdiff.^2 ./ (knots(j+1)-knots(j)); %weighted average of the square of the actual subintegral                 
                p1=ind(j);
                p2=ind(j+1);
                if p2-p1>1  %if the subinterval contains more than 1 samples
                    values=r(F,x,p1+1:p2-1,p1,p2,actFdiff); 
                    erreffective(p1+1:p2-1)=values;
                    [chosenmax(j-chosen+1),chosenind(j-chosen+1)]=max(values);
                    chosenind(j-chosen+1)=chosenind(j-chosen+1)+p1;                 
                end                
        end
        maxok(chosen+2:i+1)=maxok(chosen+1:i);
        maxok(chosen:chosen+1)=chosenmax;
        maxind(chosen+2:i+1)=maxind(chosen+1:i);
        maxind(chosen:chosen+1)=chosenind;
        wavgFdiffj(chosen+2:i+1)=wavgFdiffj(chosen+1:i);        
        wavgFdiffj(chosen:chosen+1)=chosenwavgFdiff;
        wavgFdiffall=sum(Fdiff.^2 ./ (knots(2:end)-knots(1:end-1))); %weighted average of the square of all the subintegrals 
        aprxerrj=wavgFdiffall-wavgFdiffj(1:i+1);
        
        [m,mi]=max(maxok(1:i+1)+aprxerrj);
        chosen=mi; %index of the interval, which is divided by the new knot
        maxindex = maxind(mi);
        elem=x(maxindex);
        knots=sort([knots elem]);
        ind=sort([ind maxindex]);
        %% Displaying the approximation
        if show
            [c,err,aprx]=bspline_coeffs(y,[x(1) knots x(N)] ,1);
            plot(x,y,'b',x,aprx,'r');
            legend('Original',sprintf('Approx. (PRD: %.1f%%)',err));
            drawnow;
        end
    end
    [c,aprx]=bspline_coeffs(y,[x(1) knots x(N)],1,0);
    err=norm(y-aprx)./norm(y-mean(y))*100;
end

%% Auxiliary functions

function v=r(F,x,i,p1,p2,norma)
    r1=(F(i-1)-F(p1-1)).^2./(x(i)-x(p1));
    r2=(norma-(F(i-1)-F(p1-1))).^2./(x(p2)-x(i));    
    v=r1+r2;
end
    