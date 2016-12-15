% compress_Bspline - Compress an ECG record of the MIT-BIH database. 
%
% Usage: 
%     [prd,bits,ecg,meanecg,stdecg,lens,beats,atr]=compress_Bspline(ecg_dbdir,rec_name,lead,onset,offset,outfile,order,prd_limit,acc,show);
%
% Input parameters:
%     ecg_dbdir : path of the database. 
%     recname   : recordname .
%     lead      : lead of the record.
%     onset     : begining of the signal (i.e., '00:00:00').
%     offset    : end of the signal (i.e., '00:10:00').
%     outfile   : path of the outputfile.
%     order     : order of the B-spline used for the approximations
%     prd_limit : the algorithm works till the approximation error is less than 'prd_limit'   
%     acc       : 1x2 matrix which contains the number of bits used for quantizing 
%                 the kntos and the coefficients:
%                               |  KnotsBit    CoeffBit  |
%
%     show      : optional logical value to display each step of the optimization process
%
% Output parameters:
%     prd     : approximation error of each beat in terms of PRD=(norm(sig-aprx)/norm(sig-mean(sig)))*100;
%     bits    : size of the compressed outputfile in bits.
%     ecg     : original ECG signal
%     meanecg : mean of the ecg signal
%     stdecg  : std of the ecg signal
%     lens    : length of each beat in samples 
%     beats   : cell array containingthe original, but segmented beats
%               For instance: beats{i} gives the samples of the ith beat
%     atr     : QRS annotation of the original ECG signal from 'onset' to 'offset'. 
%
% References:
% [1] M. Karczewicz, M. Gabbouj, ECG data compression by spline
%     approximation, Signal Processing, vol. 59, pp. 43-59, 1997.
%
%
%  Copyright (c) 2014, Péter Kovács <kovika@inf.elte.hu>  
%  Eötvös Lorand University, Budapest, Hungary, 2014.   
%   
%  Permission to use, copy, modify, and/or distribute this software for  
%  any purpose with or without fee is hereby granted, provided that the  
%  above copyright notice and this permission notice appear in all copies.  
%  
%  THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL  
%  WARRANTIES WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED  
%  WARRANTIES OF MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR  
%  BE LIABLE FOR ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES  
%  OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS,  
%  WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION,  
%  ARISING OUT OF OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS  
%  SOFTWARE.  

function [prd,bits,ecg,meanecg,stdecg,lens,beats,atr]=compress_Bspline(ecg_dbdir,rec_name,lead,onset,offset,outfile,order,prd_limit,acc,show)

if nargin<10
    show=0;
end

%Step 1: Segmenting the ECG signal.%
[ecg,beats,meanecg,stdecg,Rpos,atr]=getWave(rec_name,lead,ecg_dbdir,onset,offset);

%Step 2: Approximating each beat by MT-rational functions.%
knots=cell(size(beats));
coeffs=cell(size(beats));
bl_ends=zeros(size(beats),2);
l=zeros(size(beats));
bstd=zeros(size(beats));
for i=1:1:length(beats)
    strcat(num2str(i),' / ', num2str(length(beats)))
    l(i)=length(beats{i});
    [s,knots{i},coeffs{i},prd(i)]=compress(beats{i}.',order,prd_limit,1:2:length(beats{i}),show);
    bl_ends(i,:)=[beats{i}(1),beats{i}(end)];
%    [aprx prd(i)]=decompress(beats{i}',knots{i},coeff{i},order,bl_ends(i,:),1);    
end
lens=l;

%Step 3: Writting the output file.%
outfname=strcat(outfile,'.dat');
fid=fopen(outfname,'w');

%Stroing the number of beats.%
fwrite(fid,length(beats),'integer*2');

%Storing the order of Bsplines.%
fwrite(fid,order,'ubit4');

%Stroing the coefficients, knots and the first and the last values of each beat.%
knot=[];
coeff=[];
knotnum=zeros(1,length(knots));
for i=1:1:length(knots)
    knot=[knot knots{i}(2:end)-knots{i}(1:end-1)];    
    knotnum(i)=length(knots{i});
    coeff=[coeff coeffs{i}(2:end-1)'];    
end
maxcf=max(coeff);
mincf=min(coeff);
minbl=min(min(bl_ends,[],1));
maxbl=max(max(bl_ends,[],1));

%Storing statistical parameters.%
fwrite(fid,meanecg,'real*4');
fwrite(fid,stdecg,'real*4');
fwrite(fid,maxcf,'real*4');
fwrite(fid,mincf,'real*4');
fwrite(fid,maxbl,'real*4');
fwrite(fid,minbl,'real*4');

%Storing the knots for each beat.%
kn_bits=acc(1,1);
fwrite(fid,knotnum-1,strcat('ubit',num2str(kn_bits))); %knotnum-1 because we do not store the first knot, it is always equalt ot '1'.   
fwrite(fid,knot,strcat('ubit',num2str(kn_bits)));    

%Storing the coefficients for each beat.%
cf_bits=acc(1,2);
coeff_q=quant_ble(coeff,maxcf,mincf,cf_bits);
fwrite(fid,coeff_q,strcat('ubit',num2str(cf_bits)));    

%Storing the endpoints for each beat.%
ble_bits=acc(1,1);
ble=quant_ble([bl_ends(:,1); bl_ends(:,2)],maxbl,minbl,ble_bits);
fwrite(fid,ble,strcat('ubit',num2str(ble_bits)));    
fclose(fid);



% for i=1:1:length(beats)
%     %Rounding with respect to eps.
%    knot=knots{i};
%    coeff=coeffs{i}(2:end-1);
%    ble=bl_ends(i,:);
%    data=[knot,coeff',ble];
%     
%     fwrite(fid,length(data),'uint8');
%     fwrite(fid,data,'real*4');
%    fwrite(fid,0,'real*4'); %0 ends the data for each beat
% end
% fclose(fid);

%Step 3: Computing the bits that was used for the compression.%
nfo_zip=dir(strcat(outfile,'.dat'));
bits=nfo_zip.bytes*8;




%%AUXILIARY FUNCTIONS%%

%Quantizing the values of end points and coefficients.
function quantdata=quant_ble(data,maxdata,mindata,eps)
d=data-mindata;
q=(maxdata-mindata)/(2^(eps)-1);
if 0==q
    quantdata=ones(size(data));
else
    quantdata=round(d/q);
end

