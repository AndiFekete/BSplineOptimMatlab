%  Last Modified: July 25, 2014.
%  Version 1.0.
%
% decompress_Bspline - Decompress a compressed ECG record of the MIT-BIH database. 
%
% Usage: 
%     [ecg,meanecg,stdecg,knots,coeffs,bl_ends]=decompress_Bspline(fname,acc,lens)
%
% Input parameters:
%     fname : path of the compressed data. 
%     acc       : 1x2 matrix which contains the number of bits used for quantizing 
%                 the kntos and the coefficients:
%                               |  KnotsBit    CoeffBit  |
%
%     lens  : length of each beat in samples 
%
% Output parameters:
%     ecg     : reconstructed ECG signal.
%     meanecg : mean of the original ECG signal
%     stdecg  : std of the original ECG signal
%     knots   : knots of the B-splines
%     coeffs  : coefficients of the approximation
%     bl_ends : bl_ends(i,:) contains the amplitudes of the endpoints of
%               the ith beat. It is required for the reconstruction.
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

function [ecg,meanecg,stdecg,knots,coeffs,bl_ends]=decompress_Bspline(fname,acc,lens)

%Reading the input file.%
fid=fopen(strcat(fname,'.dat'),'r');
%Reading the number of beats.
beatnum=fread(fid,1,'integer*2');

%Reading the order of Bsplines.%
order=fread(fid,1,'ubit4');

knots=cell(1,beatnum);
coeffs=cell(1,beatnum);
bl_ends=zeros(beatnum,2);
l=zeros(1,beatnum);
ecg=[];

%Reading statistical parameters.%
meanecg=fread(fid,1,'real*4');
stdecg=fread(fid,1,'real*4');
maxcf=fread(fid,1,'real*4');
mincf=fread(fid,1,'real*4');
maxbl=fread(fid,1,'real*4');
minbl=fread(fid,1,'real*4');

%Reading the knots for each beat.%
kn_bits=acc(1,1);
knotnum=fread(fid,beatnum,strcat('ubit',num2str(kn_bits))); 
knot=fread(fid,sum(knotnum),strcat('ubit',num2str(kn_bits)));   

%Reading the coefficients for each beat.%
cf_bits=acc(1,2);
coeffnum=sum(knotnum)+beatnum+length(knotnum)*(order-4);
coeff=fread(fid,coeffnum,strcat('ubit',num2str(cf_bits)));
coeff_q=dequant_ble(coeff,maxcf,mincf,cf_bits);

%Reading the endpoints for each beat.%
ble_bits=acc(1,1);
ble=fread(fid,2*beatnum,strcat('ubit',num2str(ble_bits)));    
ble_q=dequant_ble(ble,maxbl,minbl,ble_bits);
bl_ends=[ble_q(1:beatnum), ble_q(beatnum+1:end)];
fclose(fid);


% eps=10^floor(log10(2^acc(1)));
% for i=1:1:beatnum
%     %Reading data for each beat.%
% %     data=[];
% %     element=fread(fid,1,'real*4');
% %     while 0~=element
% %         data=[data, element];
% %         element=fread(fid,1,'real*4');
% %     end
%     data_len=fread(fid,1,'uint8');
%     data=fread(fid,data_len,'real*4');    
%     m2=(length(data(1:end-2))-order)/2;
%     knots{i}=data(1:m2+2)';
%     coeffs{i}=[0 data(m2+3:end-2).' 0];
%     bl_ends(i,:)=data(end-1:end);
%     
    topcf=0;
    for i=1:1:beatnum  
        onseti=sum(knotnum(1:i-1))+1;
        offseti=sum(knotnum(1:i));      
        knots{i}=knot(onseti:offseti)'; %separating knots for each beat
        knots{i}=cumsum([1 knots{i}]); %restoring the knotnumber from differences 
        cnumi=length(knots{i})-4+order;
        coeffs{i}=[0 coeff_q(topcf+1:topcf+cnumi)' 0];
        topcf=topcf+cnumi;
    end

    %Reconstructing the signal.%
     for i=1:1:beatnum  
         [aprx]=decompress(rand(1,lens(i)),knots{i},coeffs{i},order,bl_ends(i,:),0);    
%         mpoles=periodize_poles(multiply_poles(p{i},ps{dbest(i)}(1:end-1)),1);
%         aprx=real(mt_generate(knots{i}(end),mpoles,c{i}));
%         slope=(bl_ends(i,2)-bl_ends(i,1))/(knots{i}(end)-1);
%         x=1:1:knots{i}(end);
%         base_line=bl_ends(i,1)+slope.*(x-1);
%         aprx=aprx+base_line;
         if 1==i
             ecg=[ecg aprx(1:end)];
         else
             ecg=[ecg aprx(2:end)];
         end
     end

 
%%AUXILIARY FUNCTIONS%%

%Reconstructing the values of end points.
function dequantdata=dequant_ble(data,maxdata,mindata,eps)
q=(maxdata-mindata)/(2^(eps)-1);
dequantdata=data*q+mindata;
