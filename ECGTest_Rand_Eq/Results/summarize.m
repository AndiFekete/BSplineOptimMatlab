%% Summarizing the results for record 117
figure(1);
hold on;
load('ResultsECGVarProUniform117.mat');
[v,indvp]=sort(prdVarPro);
plot(1:length(prdVarPro),prdVarPro(indvp),'g-','LineWidth',2);

load('ResultsECGNelderUniform117.mat');
[v,indn]=sort(prdNelder);
plot(1:length(prdNelder),prdNelder(indn),'k-','LineWidth',2);

load('ResultsECGVarProRand117.mat');
[v,indvp]=sort(prdVarPro);
plot(1:length(prdVarPro),prdVarPro(indvp),'r-','LineWidth',2);

load('ResultsECGNelderRand117.mat');
[v,indn]=sort(prdNelder);
plot(1:length(prdNelder),prdNelder(indn),'b-','LineWidth',2);

hold off;
axis([1 1515 0 80]);
box on;
h=legend('VarProUniform','NelderUniform','VarProRand','NelderRand');
set(h,'FontSize',13);
xlabel('Szívütések','FontSize',13);
ylabel('PRD (%)','FontSize',13);

%% Summarizing the results for record 119
figure(2);
hold on;
load('ResultsECGVarProUniform119.mat');
[v,indvp]=sort(prdVarPro);
plot(1:length(prdVarPro),prdVarPro(indvp),'g-','LineWidth',2);

load('ResultsECGNelderUniform119.mat');
[v,indn]=sort(prdNelder);
plot(1:length(prdNelder),prdNelder(indn),'k-','LineWidth',2);

load('ResultsECGVarProRand119.mat');
[v,indvp]=sort(prdVarPro);
plot(1:length(prdVarPro),prdVarPro(indvp),'r-','LineWidth',2);

load('ResultsECGNelderRand119.mat');
[v,indn]=sort(prdNelder);
plot(1:length(prdNelder),prdNelder(indn),'b-','LineWidth',2);


hold off;
axis([1 1965 0 90]);
box on;
h=legend('VarProUniform','NelderUniform','VarProRand','NelderRand');
set(h,'FontSize',13);
xlabel('Szívütések','FontSize',13);
ylabel('PRD (%)','FontSize',13);
