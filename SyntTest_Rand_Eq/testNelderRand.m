%% Testing fastnullad + VarPro on ECG signals

addpath('./VarPro');
addpath('./FastNullad');
addpath('./Nelder');

load('./Data/syntdata.mat','sig','knots');
load('./Results/InitRandKnots.mat','knots0');

knotsNelder=cell(size(knots));
prdnullad0=cell(size(knots));
prdnullad3=cell(size(knots));
prdNelder=cell(size(knots));

for Nkn=1:length(sig) %Nkn+5 is the number of knots.
    Nit=20; %Number of iterations of VarPro.
    options = optimset('lsqnonlin');
    options = optimset(options,'MaxIter',Nit);
    options = optimset(options,'Display','off');
    
    knotsNelder{Nkn}=zeros(size(knots{Nkn}));
    prdnullad0{Nkn}=zeros(1,size(knots{Nkn},2));
    prdnullad3{Nkn}=zeros(1,size(knots{Nkn},2));
    prdNelder{Nkn}=zeros(1,size(knots{Nkn},2));

    for i=1:1:size(knots{Nkn},2)
        display(sprintf('Processing signal pack: %d and signal: %d',length(sig),i));        
        x=(1:1:size(sig{Nkn},1))';
        y=preserve(x,sig{Nkn}(:,i));
        %Computing initial knots by fastnullad.
        [c,aprx0]=bspline_coeffs(y,knots0{Nkn}(:,i),1,false);
        prdnullad0{Nkn}(i)=norm(y-aprx0')/norm(y-mean(y))*100;
        [c,aprx]=bspline_coeffs(y,knots0{Nkn}(:,i),4,false);
        prdnullad3{Nkn}(i)=norm(y-aprx')/norm(y-mean(y))*100;
        %Improving initial knots by MatLab FminSearch.
        order=4; x=1:length(y); w=ones(size(y)); show=false; 
        initalpha=reshape(knots0{Nkn}(:,i),length(knots0{Nkn}(:,i)),1);
        initalpha=initalpha(2:end-1); %boundary knots cannot be variable.
        cf=@(alpha) costfun(order,x,y,alpha,show);
        %Optimization parameters
        alpha= fminsearch(cf,initalpha,options);
        knotsNelder{Nkn}(:,i)=[knots0{Nkn}(1,i); alpha; knots0{Nkn}(end,i)];
        [c,aprxnelder]=bspline_coeffs(y,knotsNelder{Nkn}(:,i),4,show);
        prdNelder{Nkn}(i)=norm(y-aprxnelder')/norm(y-mean(y))*100;    
    end
end
save('ResultsNelderRand.mat','knots0','knotsNelder','prdnullad0','prdnullad3','prdNelder');

