function bs=explicitBspline(x,n)
bs=zeros(size(x));
for k=0:1:n+1
    bs=bs+nchoosek(n+1,k) .* (-1).^k .* xplus( x+(n+1)/2-k ).^n; 
end
bs=bs/factorial(n);

function y=xplus(x)
for i=1:1:size(x,2)
    if x(i)>0
        y(i)=x(i);
    else
        y(i)=0;
    end 
end