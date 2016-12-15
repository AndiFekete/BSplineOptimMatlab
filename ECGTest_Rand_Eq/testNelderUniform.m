%% Testing fastnullad + VarPro on ECG signals

addpath('./VarPro');
addpath('./FastNullad');
addpath('./Nelder');

load('./Records/record117','beats');

%Parameters
Nkn=20; %Number of knots.
Nit=20; %Number of iterations of VarPro.
options = optimset;
options = optimset(options,'MaxIter',Nit);
options = optimset(options,'Display','off');

knots0=zeros(Nkn,length(beats));
knotsNelder=zeros(Nkn,length(beats));
prdnullad0=zeros(1,length(beats));
prdnullad3=zeros(1,length(beats));
prdNelder=zeros(1,length(beats));

for i=1:1:length(beats)
    display(sprintf('Processing signal %d of %d',i,length(beats))); 
    x=(1:1:length(beats{i}))';
    y=preserve(x,beats{i});
    %Computing random initial knots.
    knots0(:,i)=linspace(1,length(beats{i}),Nkn);
    [c,aprx]=bspline_coeffs(y,knots0(:,i),4,false);
    prdnullad3(i)=norm(y-aprx')/norm(y-mean(y))*100;
    %Improving initial knots by MatLab FminSearch.
    order=4; x=1:length(y); w=ones(size(y)); show=false; 
    initalpha=reshape(knots0(:,i),length(knots0(:,i)),1);
    initalpha=initalpha(2:end-1); %boundary knots cannot be variable.
    cf=@(alpha) costfun(order,x,y,alpha,show);
    %Optimization parameters
    alpha= fminsearch(cf,initalpha,options);
    knotsNelder(:,i)=[knots0(1,i); alpha; knots0(end,i)];
    [c,aprxnelder]=bspline_coeffs(y,knotsNelder(:,i),4,show);
    prdNelder(i)=norm(y-aprxnelder')/norm(y-mean(y))*100;    
end
save('ResultsECGNelderUniform117.mat','knots0','knotsNelder','prdnullad3','prdNelder');
