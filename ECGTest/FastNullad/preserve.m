function [ yy ] = preserve(x,y)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
yy=y-( (y(end)-y(1))/(x(end)-x(1)) ).*x - ((x(end)*y(1)-x(1)*y(end))/(x(end)-x(1))); 
end

