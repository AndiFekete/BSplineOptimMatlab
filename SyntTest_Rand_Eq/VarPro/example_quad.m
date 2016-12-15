%Demonstrating varpro for Bspline ECG knot optimization

%Loading data, initial knots, libraries, etc.
load('beat','beat');
load('knots','knots');

%Setting parameters.
x=linspace(0,1,100);
y=2*exp(-(x+1).^2)+3*exp(-(x-3).^2);
y=y'; x=x';
ada=@(alpha) ada_quad(alpha,x,y);
%Although we do not use weights here, it is possible to emphasize the 
%approximation of the QRS complex via weighting the data.
w=ones(size(y));
initalpha=[4; -2]; %alpha should be a column vector
n=2; %(number of knots)-1 + (degree of B-splines)
%Knots should be valid sample indices, e.g., they are integers from the interval [1,length(beat)].
lb=[-5; -5]; 
ub=[5; 5];

[alpha, c, wresid, wresid_norm, y_est, Regression] = ...
varpro(y, w, initalpha, n, ada, lb, ub);
