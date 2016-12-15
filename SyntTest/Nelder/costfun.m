function [err,aprx,co]=costfun(k,x,y,alpha,show)
if nargin<3
    show=false;
end
%% Computing the coefficients of the expansion.
sortalpha=sort(mod(alpha,length(x)-1));
knots=[x(1); sortalpha; x(end)];
[co,aprx]=bspline_coeffs(y,knots,k,show);
err=norm(y-aprx')/norm(y-mean(y))*100;

%% Displaying the results.
if show
    t=1:1:length(y);
    plot(t,y,'b',t,aprx,'r','LineWidth',2);
    hold on;    
    stem(knots,zeros(size(knots)),'k.','MarkerSize',20);
    hold off;
    h=legend('Original ECG signal',sprintf('Apprx. PRD=%.2f%%',err));
    set(h,'FontSize',14);
    axis square;
    box on;
    drawnow;
    display(sprintf('PRD=%.2f%%',err));
end