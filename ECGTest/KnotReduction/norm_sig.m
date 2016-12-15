%Normalizing the signal by baseline substraction.%
function [normsig,base_line]=norm_sig(signal)
slope=(signal(end)-signal(1))/(length(signal)-1);
x=1:1:length(signal);
base_line=signal(1)+slope.*(x-1);
normsig=signal-base_line;

