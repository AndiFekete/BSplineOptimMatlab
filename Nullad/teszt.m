%% Test data
clear all;
addpath('..');
load('beat');
x=1:length(beat);
y=beat;

%% Parameters
N=100;
timenaive=zeros(1,N);
timenullad=zeros(10,N);
timegyors=zeros(10,N);
timemin=zeros(1,N);

for i=1:10
    for k=2:1:N
        %% Naive implementation
        tic;
        [knotsnaive,err,aprx]=naive(x,y,k,0);
        timenaive(k)=toc;

        %% Effective implementation
        tic;
        [knotsnullad,err,aprx]=nullad(x,y,k,0);
        timenullad(i,k)=toc;
        
        %% Fast implementation
        tic;
        [knotsgyors,err,aprx]=gyorsnullad(x,y,k,0);
        timegyors(i,k)=toc;

%         %% Min implementation
%         tic;
%         knotsmin=nullad_min(x,y,k,0);
%         timemin(k)=toc;
        
        
%         %% Displaying the results
%         plot(1:k,timenaive(1:k),'r',1:k,timenullad(1:k),'g',1:k,timegyors(1:k),'b',1:k,timemin(1:k),'k');
%         legend('Naive','Nullad','Gyors');
%         title('Execution time of each method');    

    %    display(sprintf('l2 norm of the difference of the knots'));     
    %    display(sprintf('naive-nullad: %.2e',norm(knotsnaive-knotsnullad)));         
    %    display(sprintf('naive-gyors: %.2e',norm(knotsnaive-knotsgyors)));         

    drawnow;
    end
end