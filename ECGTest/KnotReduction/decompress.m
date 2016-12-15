function [s,prd]=decompress(beat,knots,coeff,order,bl,show)
x=1:1:length(beat);
s=zeros(1,length(beat));
tt=[knots(1).*ones(1,order-1) knots knots(end).*ones(1,order-1)];

%Bspline reconstruction.%
for i=1:1:length(coeff)
    p=coeff(i)*Bspline(order-1,i,tt,x);
    s=s+p;
%     plot(x,c(i)*Bspline(order-1,i,tt,x))
%     drawnow
%     hold on
%     pause
end
%Baseline reconstruction.%
slope=(bl(end)-bl(1))/(length(beat)-1);
base_line=bl(1)+slope.*(x-1);
s=s+base_line;

prd=sqrt(sum((beat-s).^2)/sum((beat-mean(beat)).^2))*100;
cr=(length(knots)+length(coeff)-1)/length(beat)*100;

if show
    plot(x,beat,'b',x,s,'r','LineWidth',2);
    legend(sprintf('CR: %.1f',cr),sprintf('PRD: %.1f',prd));
end