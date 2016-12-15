function [Phi,dPhi,Ind] = ada_bspline(k,x,y,alpha,show)
%
%     This function computes Phi and, if possible, dPhi.
%
%     On Input: 
%
%        alpha q x 1    contains the current value of the alpha parameters.
%
%        Note:  If more input arguments are needed, use the standard
%               Matlab syntax to accomplish this.  For example, if
%               the input arguments to ada are t, z, and alpha, then
%               before calling varpro, initialize t and z, and in calling 
%               varpro, replace "@ada" by "@(alpha)ada(t,z,alpha)".
%
%     On Output:
%
%        Phi   m x n1   where Phi(i,j) = phi_j(alpha,t_i).
%                       (n1 = n if there is no extra term; 
%                        n1 = n+1 if an extra term is used)
%        dPhi  m x p    where the columns contain partial derivative
%                       information for Phi and p is the number of 
%                       columns in Ind 
%                       (or dPhi = [] if derivatives are not available).
%        Ind   2 x p    Column k of dPhi contains the partial
%                       derivative of Phi_j with respect to alpha_i, 
%                       evaluated at the current value of alpha, 
%                       where j = Ind(1,k) and i = Ind(2,k).
%                       Columns of dPhi that are always zero, independent
%                       of alpha, need not be stored. 
%        Example:  if  phi_1 is a function of alpha_2 and alpha_3, 
%                  and phi_2 is a function of alpha_1 and alpha_2, then 
%                  we can set
%                          Ind = [ 1 1 2 2
%                                  2 3 1 2 ]
%                  In this case, the p=4 columns of dPhi contain
%                          d phi_1 / d alpha_2,
%                          d phi_1 / d alpha_3,
%                          d phi_2 / d alpha_1,
%                          d phi_2 / d alpha_2,
%                  evaluated at each t_i.
%                  There are no restrictions on how the columns of
%                  dPhi are ordered, as long as Ind correctly specifies
%                  the ordering.
%
%        If derivatives dPhi are not available, then set dPhi = Ind = [].

[sortalpha,sortind]=sort(alpha);
Phi=calcmat(k,x,[x(1); sortalpha; x(end)]);
[dPhi, Ind]=dercalcmat(k,x,[x(1); sortalpha; x(end)],sortind);
% dPhi_aprx=zeros(size(dPhi));
% h=0.001;
% for i=1:1:size(dPhi,2)
%     Phi1=Phi(:,Ind(1,i));
%     xk=alpha;
%     xk(Ind(2,i))=xk(Ind(2,i))+h;
%     Phi2=BsplineDeBoor(k-1,Ind(1,i)+1,[x(1)*ones(k,1); sort(xk); x(end)*ones(k,1)],x)';
%     dPhi_aprx(:,i)=(Phi2-Phi1)/h;
% end
%dPhi=[]; Ind=[];

%Displaying the results at each step.
if show
    [c,aprx]=bspline_coeffs(y,[x(1); alpha; x(end)],k,false);
    knotvals=interp1(x,aprx,[x(1); alpha; x(end)],'spline',0);
    prd=norm(y-aprx')/norm(y-mean(y))*100;
    display(sprintf('PRD of the VarPro approximation: %.2f%%',prd));
    plot(x,y,'b',x,aprx,'r','LineWidth',2);
    hold on;
    stem([x(1); alpha; x(end)],knotvals,'r');
    hold off;
    h=legend('ECG',sprintf('VarPro PRD: %.2f%%',prd));
    set(h,'FontSize',13);
    grid on;
    axis tight;
    drawnow;
    pause(0.5);
end