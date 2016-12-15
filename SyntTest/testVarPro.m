%% Testing fastnullad + VarPro on ECG signals

addpath('./VarPro');
addpath('./FastNullad');

load('./Data/syntdata.mat','sig','knots');

knots0=cell(size(knots));
knotsVarPro=cell(size(knots));
prdnullad0=cell(size(knots));
prdnullad3=cell(size(knots));
prdVarPro=cell(size(knots));

for Nkn=1:length(sig) %Nkn+5 is the number of knots.
    Nit=20; %Number of iterations of VarPro.
    options = optimset('lsqnonlin');
    options = optimset(options,'MaxIter',Nit);

    knots0{Nkn}=zeros(size(knots{Nkn}));
    knotsVarPro{Nkn}=zeros(size(knots{Nkn}));
    prdnullad0{Nkn}=zeros(1,size(knots{Nkn},2));
    prdnullad3{Nkn}=zeros(1,size(knots{Nkn},2));
    prdVarPro{Nkn}=zeros(Nkn,size(knots{Nkn},2));

    for i=1:1:size(knots{Nkn},2)
        display(sprintf('Processing signal pack: %d and signal: %d',length(sig),i));        
        x=(1:1:size(sig{Nkn},1))';
        y=preserve(x,sig{Nkn}(:,i));
        %Computing initial knots by fastnullad.
        [knots0{Nkn}(:,i),aprx,err]=gyorsnullad(x',y',Nkn+5-2,false);
        [c,aprx0]=bspline_coeffs(y,knots0{Nkn}(:,i),1,false);
        prdnullad0{Nkn}(i)=norm(y-aprx0')/norm(y-mean(y))*100;
        [c,aprx]=bspline_coeffs(y,knots0{Nkn}(:,i),4,false);
        prdnullad3{Nkn}(i)=norm(y-aprx')/norm(y-mean(y))*100;
        %Improving initial knots by VarPro.
        order=4; x=1:length(y); w=ones(size(y)); show=false; 
        ada=@(alpha) ada_bspline(order,x,y,alpha,show);    
        initalpha=reshape(knots0{Nkn}(:,i),length(knots0{Nkn}(:,i)),1);
        initalpha=initalpha(2:end-1); %boundary knots cannot be variable.
        n=length(initalpha)-1 + (order-1); %(number of knots)-1 + (degree of B-splines)
        %Knots should be valid sample indices, e.g., they are integers from the interval [1,length(beat)].
        lb=ones(length(initalpha),1)+1; 
        ub=size(sig{Nkn},1)*ones(length(initalpha),1)-1;
        [alpha, c, wresid, wresid_norm, y_est, Regression] = ...
        varpro(y, w, initalpha, n, ada, lb, ub, options);
        knotsVarPro{Nkn}(:,i)=[knots0{Nkn}(1,i); alpha; knots0{Nkn}(end,i)];
        [c,aprxvarpro]=bspline_coeffs(y,knotsVarPro{Nkn}(:,i),4,false);
        prdVarPro{Nkn}(i)=norm(y-aprxvarpro')/norm(y-mean(y))*100;    
    end
end
save('ResultsVarPro.mat','knots0','knotsVarPro','prdnullad0','prdnullad3','prdVarPro');
