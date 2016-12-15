%Generating sygnals for synthetic B-spline tests
N=360; %Number of samples
n=25;  %Maximum number of knots
M=100; %Number of generated signals for each knot configuration
a=-20; b=20; %Interval for the coefficients.
order=4;     %Order of the B-splines

knots=cell(1,20);
coeffs=cell(1,20);
sig=cell(1,20);
x=1:1:N;
for ni=6:n
    knots{ni-5}=zeros(ni,M);
    coeffs{ni-5}=zeros(ni,M);
    sig{ni-5}=zeros(N,M);
    for mi=1:1:M        
        innerknots=sort(randperm(N-2,ni-2)+1);            
        knots{ni-5}(:,mi)=[1; innerknots'; N];
        coeffs{ni-5}(:,mi)=a+rand(ni,1)*(b-a);
        sig{ni-5}(:,mi)=generateBspline(x,coeffs{ni-5}(:,mi),knots{ni-5}(:,mi),order);
        %plot(x,sig{ni-5}(:,mi),'r',knots{ni-5}(:,mi),sig{ni-5}(knots{ni-5}(:,mi),mi),'b*');
    end
end
save('syntdata.mat','knots','coeffs','sig');


