%--------------------Example of the preserved signal----------------------%
function [s,knots,coeff,prd]=compress(sig,order,knot_limit,init_knot,show)
x=1:length(sig);
psig=preserve(x,sig);

knots=init_knot;
all_knots=knots;
not_knots=[];
[c,prd,s,Rl,zl]=bspline_coeffs(psig,knots,order);

step=1;
%saveas(gca, strcat('./anim2/bspline_ecg_',int2str(step),'.jpg'));

%Removing knots.%
while length(knots)>knot_limit % prd<prd_limit
    wpp=zeros(1,length(knots)-2);
    for i=2:length(knots)-1
        wpp(i-1)=predict_mse(psig,s,c,order,knots,i); 
    end
    [vv,indd]=min(wpp);
    not_knots(end+1)=knots(indd+1);
    tt=knots;
    knots(indd+1)=[];
    [c,prd,s,Rl,zl]=bspline_coeffs(psig,knots,order,indd+1,tt,Rl,zl,show);
    if show
        hold on;
        width=(max(psig)-min(psig));
        bl=min(s)-0.1*width;
        plot(not_knots,bl.*ones(1,length(not_knots)),'bx','MarkerSize',12,'LineWidth',2);
        plot(knots,bl.*ones(1,length(knots)),'r.','MarkerSize',15);
        legend(sprintf('PRD: %.02f%%',prd),sprintf('CR: %.0f%%',100*(length(c) + length(knots))/(length(sig))));
        axis([0 length(sig),min(psig)-0.2*width,max(psig)+0.1*width]);
        hold off;
        drawnow;
        step=step+1;
        %saveas(gca, strcat('./anim2/bspline_ecg_',int2str(step),'.jpg'));        
    end
end
coeff=c;
%------------------------Plotting final results.--------------------------%
if show
    subplot(2,2,[1 4])
    [c,prd,s]=bspline_coeffs(psig,knots,order);
    plot(s,'r','LineWidth',2);
    hold on;
    plot(psig,'g');    
    width=(max(s)-min(s));
    bl=min(s)-0.1*width;
    plot(not_knots,bl.*ones(1,length(not_knots)),'bx','MarkerSize',12,'LineWidth',2);
    stem(knots,psig(knots),'r.','MarkerSize',15);
    legend(sprintf('PRD: %.02f%%',prd));
    axis([0 length(sig),min(psig)-0.2*width,max(psig)+0.1*width]);
    title('Approximation');
    hold off;
end