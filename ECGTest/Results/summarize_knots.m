%% Displaying the knot vectors
%---------------------------fastnullad knots------------------------------%
load('ResultsECGNelder119.mat');
load('ResultsECGKnotRed119.mat');

load('ResultsECGVarPro119.mat');
knotnum=20;
figure(1); 
%subplot(1,3,1);
hold on;
for i=1:150%size(knots0,2) %the first two elemnts of matdir are '.' and '..'
    plot(knots0(:,i),i,'r.');
end
hold off;
box on;
%title('Nullad knots');
xlabel('Alappontok indexe','FontSize',13);
ylabel('Szívütések indexe','FontSize',13);
axis tight;

%---------------------------VarPro knots------------------------------%
load('ResultsECGVarPro119.mat');
knotnum=20;
figure(2); 
%subplot(1,3,2);
hold on;
for i=1:150%size(knotsVarPro,2) %the first two elemnts of matdir are '.' and '..'
    plot(knotsVarPro(:,i),i,'g.');
end
hold off;
box on;
%title('VarPro knots');
xlabel('Alappontok indexe','FontSize',13);
ylabel('Szívütések indexe','FontSize',13);
axis tight;

%---------------------------KnotRed knots------------------------------%
load('ResultsECGKnotRed119.mat');
knotnum=20;
figure(3); 
%subplot(1,3,3);
hold on;
for i=1:150%size(knotsRed,2) %the first two elemnts of matdir are '.' and '..'
    plot(knotsRed(:,i),i,'m.');
end
hold off;
box on;
%title('KnotRed knots');
xlabel('Alappontok indexe','FontSize',13);
ylabel('Szívütések indexe','FontSize',13);
axis tight;