%min keresés
function knots=nullad_min(x,y,k,show)
    N=length(x);
    knots=[2,N-1];
    n=length(knots)
    
    %k darab alappontot hozzáadunk
    for j=1:k
       error=Inf*ones(1,N); %error init
       %error kiszámolása a bels? pontokban
       for i=3:N-2        
           if (~(any(i==knots)))
                %ha már alappont, nem vizsgáljuk
                tmp=sort([knots i]);
                A=zeros(N,n);
                d=ones(1,n);
                for ii=1:n   
                    d(ii)=1/(tmp(ii+1)-tmp(ii));
                    for jj=1:N
                        A(jj,ii)=(tmp(ii)<=jj && jj<tmp(ii+1));
                    end
                end
                G=diag(d);
                c=G*A'*y';
                %error kiszámolása
                error(i)=sum((y'-A*c).^2);
           end
       end
       %
       [m,index]=min(error);
       knots=sort([knots index]);
       n=n+1;
    end
    
    if show
        [c,err,aprx]=bspline_coeffs(y,[x(1) knots x(N)] ,1);
        plot(x,y,'b',x,aprx,'r');
        legend('Original',sprintf('Approx. (PRD: %.1f%%)',err));
        drawnow;
    end
end