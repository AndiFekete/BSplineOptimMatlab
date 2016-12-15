load('ResultsVarProUniform.mat');
load('ResultsNelderUniform.mat');
load('../Data/syntdata','knots');


meanprdVarPro=zeros(1,20);
meanprdNelder=zeros(1,20);
meanprdVarProRand=zeros(1,20);
meanprdNelderRand=zeros(1,20);

meanknoterrVarPro=zeros(1,20);
meanknoterrNelder=zeros(1,20);
meanknoterrVarProRand=zeros(1,20);
meanknoterrNelderRand=zeros(1,20);

for i=1:1:20
    meanprdVarPro(i)=mean(prdVarPro{i});
    meanprdNelder(i)=mean(prdNelder{i});
    
    meanknoterrVarPro(i)=mean(mean(abs(knots{i}(2:end-1,:)-knotsVarPro{i}(2:end-1,:))));
    meanknoterrNelder(i)=mean(mean(abs(knots{i}(2:end-1,:)-knotsNelder{i}(2:end-1,:))));
end

load('ResultsVarProRand.mat');
load('ResultsNelderRand.mat');

for i=1:1:20
    meanprdVarProRand(i)=mean(prdVarPro{i});
    meanprdNelderRand(i)=mean(prdNelder{i});
    
    meanknoterrVarProRand(i)=mean(mean(abs(knots{i}(2:end-1,:)-knotsVarPro{i}(2:end-1,:))));
    meanknoterrNelderRand(i)=mean(mean(abs(knots{i}(2:end-1,:)-knotsNelder{i}(2:end-1,:))));
end

figure(1);
hold on;
plot(6:25,meanprdVarPro,'g',6:25,meanprdNelder,'k','LineWidth',2);
plot(6:25,meanprdVarProRand,'r',6:25,meanprdNelderRand,'b','LineWidth',2);
hold off;
box on;
h=legend('VarProUniform','NelderUniform','VarProRand','NelderRand');
set(h,'FontSize',13);
xlabel('Alappontok száma','FontSize',13);
ylabel('Átlagos PRD (%)','FontSize',13);
axis([5 26 0 65]);

figure(2);
hold on;
plot(6:25,meanknoterrVarPro,'g',6:25,meanknoterrNelder,'k','LineWidth',2);
plot(6:25,meanknoterrVarProRand,'r',6:25,meanknoterrNelderRand,'b','LineWidth',2);
hold off;
box on;
h=legend('VarProUniform','NelderUniform','VarProRand','NelderRand');
set(h,'FontSize',13);
xlabel('Alappontok száma','FontSize',13);
ylabel('Átlagos eltérés','FontSize',13);
axis([5 26 20 65]);
