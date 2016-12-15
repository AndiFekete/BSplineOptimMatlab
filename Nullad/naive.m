function [knots, aprx, err] = naive(x,y,k,show)
    N=length(y);
    knots=[2,N-1];
    for j=1:1:k
        err=Inf*ones(1,N);
            for i=3:N-2
               if ~sum(i==knots)
                   t=sort([knots, i]);
                   [c,err(i),aprx]=bspline_coeffs(y,[1 t N],1);        
               end
            end
        [v,ind]=min(err);
        knots=sort([knots, ind]);
        if show
            [c,e,aprx]=bspline_coeffs(y,[1 knots N],1);  
            subplot(1,2,1);
            plot(x,y,'g',x,aprx,'r');   
            subplot(1,2,2);
            plot(err(err~=0));   
        end
    end
    [c,err,aprx]=bspline_coeffs(y,[1 knots N],1); 
end

