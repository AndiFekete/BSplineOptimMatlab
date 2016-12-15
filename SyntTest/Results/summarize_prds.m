load('ResultsVarPro');
load('ResultsNelder');
load('ResultsKnotRed');
load('../Data/syntdata','knots');

meanprdnullad0=zeros(1,20);
meanprdnullad3=zeros(1,20);
meanprdVarPro=zeros(1,20);
meanprdNelder=zeros(1,20);
meanprdKnotRed=zeros(1,20);

meanknoterrnullad0=zeros(1,20);
meanknoterrVarPro=zeros(1,20);
meanknoterrNelder=zeros(1,20);
meanknoterrKnotRed=zeros(1,20);

for i=1:1:20
    meanprdnullad0(i)=mean(prdnullad0{i});
    meanprdnullad3(i)=mean(prdnullad3{i});
    meanprdVarPro(i)=mean(prdVarPro{i});
    meanprdNelder(i)=mean(prdNelder{i});
    meanprdKnotRed(i)=mean(prdKnotRed{i});
    
    meanknoterrnullad0(i)=mean(mean(abs(knots{i}(2:end-1,:)-knots0{i}(2:end-1,:))));
    meanknoterrVarPro(i)=mean(mean(abs(knots{i}(2:end-1,:)-knotsVarPro{i}(2:end-1,:))));
    meanknoterrNelder(i)=mean(mean(abs(knots{i}(2:end-1,:)-knotsNelder{i}(2:end-1,:))));
    meanknoterrKnotRed(i)=mean(mean(abs(knots{i}(2:end-1,:)-knotsRed{i}(2:end-1,:))));
end
figure(1);
hold on;
plot(6:25,meanprdnullad0,'r',6:25,meanprdnullad3,'b','LineWidth',2);
plot(6:25,meanprdVarPro,'g',6:25,meanprdNelder,'k','LineWidth',2);
plot(6:25,meanprdKnotRed,'m','LineWidth',2);
hold off;
box on;
h=legend('Nulad0','Nullad3','VarPro','Nelder','KnotRed');
set(h,'FontSize',13);
xlabel('Alappontok száma','FontSize',13);
ylabel('Átlagos PRD (%)','FontSize',13);
axis([5 26 0 80]);

figure(2);
hold on;
plot(6:25,meanknoterrnullad0,'r','LineWidth',2);
plot(6:25,meanknoterrVarPro,'g',6:25,meanknoterrNelder,'k','LineWidth',2);
plot(6:25,meanknoterrKnotRed,'m','LineWidth',2);
hold off;
box on;
h=legend('Nulad0','VarPro','Nelder','KnotRed');
set(h,'FontSize',13);
xlabel('Alappontok száma','FontSize',13);
ylabel('Átlagos eltérés','FontSize',13);
%axis([5 26 0 80]);
