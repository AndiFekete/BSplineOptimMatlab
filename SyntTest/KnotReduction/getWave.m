% getWave - Segmenting ECG signals into heart beats.%
%
% Usage: 
%     [ecg,beats,Rpos] = getWave(rec_name,fs,lead,ecg_dbdir)
%
% Input parameters:
%     rec_name  : record name related to Physionet/MIT-BIH Arrhythmia database  
%     lead      : specify the channel number to be segmented
%     ecg_dbdir : path of the records
%     onset     : starting point of the segmentation (default: '00:00:00')
%     offset    : ending point of the segmentation (default: '00:05:00')
%
% Output parameters:
%     ecg     : extracted ECG signal  
%     beats   : extracted beats of the signal 
%     meanecg : mean of the ECG signal
%     stdecg  : std of the ECG signal 
%     Rpos    : positions of the R peaks
%
%NOTE: WFDB PhysioToolkit package should be loaded in MATLAB!!!
%
%  Copyright (c) 2014, Péter Kovács <kovika@inf.elte.hu>  
%  Eötvös Lorand University, Budapest, Hungary, 2014.   

function [ecg,beats,meanecg,stdecg,Rpos,atr]=getWave(rec_name,lead,ecg_dbdir,onset,offset)

%Initializing default parameters.
if nargin<4
    onset='00:00:00';
end
if nargin<5
    offset='00:05:00';
end

% %Loading signal directory. (with deprecated wfdb library)%
% wp=what;
% cd(ecg_dbdir);
% ecg=rdsamp(rec_name,'begin',onset,'stop',offset,'phys',true);
% ecg=ecg(:,lead+1);

%Loading signal directory. (with actual version of wfdb library)%
wp=what;
cd(ecg_dbdir);
FS=360; %Sampling frequency of MIT-BIH records
N=FS*(str2num(offset(1:2))*3600+str2num(offset(4:5))*60+str2num(offset(7:8)));
N0=FS*(str2num(onset(1:2))*3600+str2num(onset(4:5))*60+str2num(onset(7:8)));
[tm,ecg,Fs]=rdsamp(rec_name,1,N,N0,1);


%Normalizing the signal.%
meanecg=mean(ecg);
ecg=ecg-meanecg;
stdecg=std(ecg);
ecg=ecg/stdecg;

% %Loading signal annotations with deprecated wfdb library.%
% atr=rdann(rec_name,strcat('atr',num2str(lead)),'start',onset,'stop',offset);
% markers_num=[atr.sampleNumber];
% if markers_num(1)>length(ecg)
%     markers_num=markers_num-length(ecg);
% end
% markers_type=[atr.typeMnemonic];
% beat_annot=['N','L','R','B','A','a','J','S','V','r','F','e','j','n','E','/','f','Q','?'];

%Loading signal annotations with the actual version of wfdb library.%
[ann,type,subtype,chan,num,comments]=rdann(rec_name,strcat('atr',num2str(lead)),[],N,N0);
markers_num=ann-N0;
markers_type=type;
atr=[]; atr.ann=markers_num; atr.type=markers_type; 
atr.subtype=subtype; atr.chan=chan; atr.num=num;
beat_annot=['N','L','R','B','A','a','J','S','V','r','F','e','j','n','E','/','f','Q','?'];

%Localizing the positions of R peaks.%
[wmark_R atr]=find_markers(atr,beat_annot);
Rpos=markers_num(wmark_R);

%Segmenting the ECG into beats.%
cut=130;
beats=cell(ceil(length(Rpos))-2,1);
for i=2:1:length(Rpos)-1
    on=max([0,Rpos(i)-cut]);
    off=Rpos(i+1)-cut;
    beats{i-1}=ecg(on:off);
end
on=max([0,Rpos(2)-cut]);
off=Rpos(length(Rpos))-cut;
ecg=ecg(on:end);
ecg=ecg(1:off-on+1);
%Deleting the annotations of the first and second beats.
atr.ann(1)=[];
atr.type(1)=[];
atr.subtype(1)=[];
atr.chan(1)=[];
atr.num(1)=[];
atr.ann(end)=[];
atr.type(end)=[];
atr.subtype(end)=[];
atr.chan(end)=[];
atr.num(end)=[];
%Correcting the annotations with the skipped first beat.
atr.ann=atr.ann-on+1;
%beats{end}=ecg(Rpos(i+1)-cut:end);
%wrann(atr, rec_name, 'wdd')
wrann(rec_name,'wdd',atr.ann,atr.type,atr.subtype,atr.chan,atr.num);
%Loading the actual working directory.%
cd(wp.path);



%%AUXILIARY FUNCTIONS.%%

%Selecting the specified wave markers(p N t) in 'markers_type'.%
function [wmark newatr]=find_markers(atr,wave_name)
top=0;
topind=0;
newatr=atr;
markers_type=[atr.type];
wmark=zeros(1,length(markers_type));
ind=ones(1,length(markers_type));
for i=1:1:length(markers_type)
    if isin(markers_type(i),wave_name)        
        top=top+1;       
        wmark(top)=i;
    else
        topind=topind+1;
        ind(topind)=i;
    end
end
wmark=wmark(1:top);
ind=ind(1:topind);
newatr.ann(ind)=[];
newatr.type(ind)=[];
newatr.subtype(ind)=[];
newatr.chan(ind)=[];
newatr.num(ind)=[];

%Counting beats in the ECG signal.%
function beat_number=count_beats(markers_type)
beat_number=0;
for i=1:1:length(markers_type)
    if markers_type(i)=='N'
        beat_number=beat_number+1;
    end
end

%Matching the marker symbol 'a' with the members of the 'set'.%
function l=isin(a,set)
l=false;
for i=1:1:length(set)
    if strcmp(a,set(i))
        l=true;
        break;
    end
end

