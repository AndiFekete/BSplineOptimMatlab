function [Phi,dPhi,Ind] = ada_bspline1knot(alpha,x,y)
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

Phi=zeros(length(x),2);
Phi(:,1)=exp(-(x+alpha(1)).^2);
Phi(:,2)=exp(-(x+alpha(2)).^2);

dPhi=zeros(length(x),2);

dPhi(:,1)=-exp(-(x+alpha(1)).^2).*(x+alpha(1))*2;
dPhi(:,2)=-exp(-(x+alpha(2)).^2).*(x+alpha(2))*2;
Ind=[1 2;1 2];
%dPhi=[]; Ind=[];
c=Phi\y;
aprx=Phi*c;
plot(x,y,'b',x,aprx,'r');
drawnow;
pause(0.5)
