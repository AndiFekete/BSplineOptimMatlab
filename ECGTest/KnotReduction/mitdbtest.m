%% Loading libraries.
addpath('C:\Users\kovika\Documents\MATLAB\CR_TEST\WFDB');
clear;

load('results','test_records','cr_basic','prd_basic');

%% Setting parameters.
for ind=[13, 14];
    ecg_dbdir='./mitdb';
    rec_name=test_records(ind,:);
    prd_limit=prd_basic(ind)-0.5;
    order=3;
    acc=[8,8];
    lead=1;
    onset='00:00:00';
    offset='00:30:00';
    outfile=['./mitdbtest/Bspline/DAT/',rec_name,'_bspline_',num2str(order)];
    show=false;

    % De/Compressing the ECG signal by using the basic algorithm.
    [prd,bits,ecgsig,meanecg,stdecg,lens,beats,atr]=compress_Bspline(ecg_dbdir,rec_name,lead,onset,offset,outfile,order,prd_limit,acc,show);
    [aprx,meanecg,stdecg]=decompress_Bspline(outfile,acc,lens);
    save(['./mitdbtest/Bspline/MAT/',rec_name,'_bspline_',num2str(order)],'prd','bits','ecgsig','meanecg','stdecg','lens','beats','aprx','atr');
    tm=0:1:length(ecgsig)-1;
    wrsamp(tm',ecgsig,[rec_name,'_bspline'],360,200,'212');
    wrann([rec_name,'_bspline'],'atr',atr.ann,atr.type,atr.subtype,atr.chan,atr.num);
    %[ann,type,subtype,chan,num]=rdann([rec_name,'_bspline'],'atr');
end