function [sig]=generateBspline(x,co,knots,order)
    tt=[knots(1).*ones(order-1,1); knots; knots(end).*ones(order-1,1)]; 
    c=[0; co; 0]; 
    sig=zeros(size(x));
    for i=1:1:length(c)
        p=c(i)*BsplineDeBoor(order-1,i,tt,x);
        sig=sig+p;
    end   
end