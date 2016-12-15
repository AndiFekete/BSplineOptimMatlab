%% Testing Knot reduction on ECG signals

addpath('./KnotReduction');

load('./Records/record117','beats');

%Parameters
Nkn=20; %Number of knots.

knotsRed=zeros(Nkn,length(beats));
prdKnotRed=zeros(1,length(beats));
for i=1:1:length(beats)
    display(sprintf('Processing signal %d of %d',i,length(beats)));    
    x=(1:1:length(beats{i}))';
    y=preserve(x,beats{i});
    init_knot=[1:2:length(beats{i})-1 length(beats{i})];
    [aprx,knotsRed(:,i),coeff,prd]=compress(y',3,Nkn,init_knot,false);
    prdKnotRed(i)=norm(y-aprx')/norm(y-mean(y))*100;    
end
save('ResultsECGKnotRed.mat','knotsRed','prdKnotRed');