%Demonstrating varpro for Bspline ECG knot optimization

%Loading data, initial knots, libraries, etc.
load('beat','beat');
load('knots','knots');

%Setting parameters.
y=reshape(beat,length(beat),1);
x=1:length(beat);
tt=sort([x(1); x(50); x(140); x(200); x(end)]);
y=3*BsplineDeBoor(3,1,tt,x)';
%ada function for evaluating basis functions and their partial derivatives.
order=4; x=1:length(y); 
ada=@(alpha) ada_bspline1knot(order,x,y,alpha);
%Although we do not use weights here, it is possible to emphasize the 
%approximation of the QRS complex via weighting the data.
w=ones(size(y));
initalpha=10; %alpha should be a column vector
n=1; %(number of knots)-1 + (degree of B-splines)
%Knots should be valid sample indices, e.g., they are integers from the interval [1,length(beat)].
lb=ones(length(initalpha),1)+1; 
ub=length(beat)*ones(length(initalpha),1)-1;

[alpha, c, wresid, wresid_norm, y_est, Regression] = ...
varpro(y, w, initalpha, n, ada, lb, ub);

%B-spline approximation
[init_c,init_aprx]=bspline_coeffs(y,[x(1); knots'; x(end)],order,false);
[c,aprx]=bspline_coeffs(y,[x(1); alpha; x(end)],order,false);
[round_c,round_aprx]=bspline_coeffs(y,round([x(1); alpha; x(end)]),order,false);
prd_init=norm(beat-init_aprx)/norm(beat-mean(beat))*100;
prd=norm(beat-aprx)/norm(beat-mean(beat))*100;
prd_round=norm(beat-round_aprx)/norm(beat-mean(beat))*100;
display(sprintf('PRD of the initial approximation: %.2f%%',prd_init));
display(sprintf('PRD of the VarPro approximation: %.2f%%',prd));
display(sprintf('PRD of the Rounded VP approximation: %.2f%%',prd_round));

plot(x,y,'b',x,init_aprx,'r',x,aprx,'g',x,round_aprx,'m','LineWidth',2);
legend('ECG',sprintf('Init PRD: %.2f%%',prd_init),sprintf('VarPro PRD: %.2f%%',prd),sprintf('Rounded VP PRD: %.2f%%',prd_round));
