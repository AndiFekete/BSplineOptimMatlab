%% Summarizing the results for record 117
figure(1);
hold on;
load('ResultsECGVarPro117.mat');
[v,ind0]=sort(prdnullad0);
[v,ind3]=sort(prdnullad3);
[v,indvp]=sort(prdVarPro);
plot(1:length(prdnullad0),prdnullad0(ind0),'r-','LineWidth',2);
plot(1:length(prdnullad3),prdnullad3(ind3),'b-','LineWidth',2);
plot(1:length(prdVarPro),prdVarPro(indvp),'g-','LineWidth',2);

load('ResultsECGNelder117.mat');
[v,indn]=sort(prdNelder);
plot(1:length(prdNelder),prdNelder(indn),'k-','LineWidth',2);

load('ResultsECGKnotRed117.mat');
[v,indkn]=sort(prdKnotRed);
plot(1:length(prdKnotRed),prdKnotRed(indkn),'m-','LineWidth',2);
hold off;
axis([1 1515 0 35]);
box on;
h=legend('Nulad0','Nullad3','VarPro','Nelder','KnotRed');
set(h,'FontSize',13);
xlabel('Szívütések','FontSize',13);
ylabel('PRD (%)','FontSize',13);

%% Summarizing the results for record 119
figure(2);
hold on;
load('ResultsECGVarPro119.mat');
[v,ind0]=sort(prdnullad0);
[v,ind3]=sort(prdnullad3);
[v,indvp]=sort(prdVarPro);
plot(1:length(prdnullad0),prdnullad0(ind0),'r-','LineWidth',2);
plot(1:length(prdnullad3),prdnullad3(ind3),'b-','LineWidth',2);
plot(1:length(prdVarPro),prdVarPro(indvp),'g-','LineWidth',2);

load('ResultsECGNelder119.mat');
[v,indn]=sort(prdNelder);
plot(1:length(prdNelder),prdNelder(indn),'k-','LineWidth',2);

load('ResultsECGKnotRed119.mat');
[v,indkn]=sort(prdKnotRed);
plot(1:length(prdKnotRed),prdKnotRed(indkn),'m-','LineWidth',2);
hold off;
axis([1 1965 0 35]);
box on;
h=legend('Nulad0','Nullad3','VarPro','Nelder','KnotRed');
set(h,'FontSize',13);
xlabel('Szívütések','FontSize',13);
ylabel('PRD (%)','FontSize',13);
