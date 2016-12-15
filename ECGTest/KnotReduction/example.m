%% Loading libraries.
addpath('D:\Kovacs_Peter\Matlab\ECG_Comp_cikk\Tests\WFDB');

%% Setting parameters.
ecg_dbdir='./records';
rec_name='117';
order=3;
prd_limit=8.0;
acc=[8,8];
lead=1;
onset='00:00:00';
offset='00:00:10';
outfile=['./records/compressed/',rec_name,'_comp'];
show=false;

%% Compressing the ECG signal.
[prd,bits,ecgsig,meanecg,stdecg,lens]=compress_Bspline(ecg_dbdir,rec_name,lead,onset,offset,outfile,order,prd_limit,acc,show);
[aprx,meanecg,stdecg,knots,coeffs,bl_ends]=decompress_Bspline(outfile,acc,lens);

%% Saving one channel of the ECG to compute the CR.
tm=0:1:length(ecgsig)-1;
wrsamp(tm',ecgsig,['./records/compressed/',rec_name,'_uncomp'],360,200,'212');
nfo_dat=dir(['./records/compressed/',rec_name,'_uncomp.dat']);
uncompbits=nfo_dat.bytes*8;
nfo_dat=dir([outfile,'.dat']);
compbits=nfo_dat.bytes*8;

%% Displaying statistical informations.
disp(sprintf('Mean PRD: %.2f %%',mean(prd)));
disp(sprintf('Compression ratio: 1:%.2f',uncompbits/compbits));

%% Plotting the original ECG and the reconstructed signal.
tm=[0:1:length(ecgsig)-1]/360;
plot(tm/360,ecgsig,'b',tm/360,aprx,'r');
legend('Original ECG','Approximation');
