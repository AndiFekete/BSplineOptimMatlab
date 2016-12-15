%% Testing Knot reduction on ECG signals

addpath('./KnotReduction');

load('./Data/syntdata.mat','sig','knots');

knotsRed=cell(size(knots));

%Parameters
Nkn=20; %Number of knots.
prdKnotRed=cell(size(knots));

for Nkn=1:length(sig) %Nkn+5 is the number of knots.
    knotsRed{Nkn}=zeros(size(knots{Nkn}));
    for i=1:1:size(knots{Nkn},2)        
        display(sprintf('Processing signal pack: %d and signal: %d',length(sig),i));
        x=(1:1:size(sig{Nkn},1))';
        y=preserve(x,sig{Nkn}(:,i));        
        init_knot=[1:2:size(sig{Nkn},1)-1 size(sig{Nkn},1)];
        [aprx,knotsRed{Nkn}(:,i),coeff,prd]=compress(y',3,Nkn+5,init_knot,false);
        prdKnotRed{Nkn}(i)=norm(y-aprx')/norm(y-mean(y))*100;    
    end
end