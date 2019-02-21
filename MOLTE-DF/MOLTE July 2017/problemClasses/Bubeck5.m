function [f, beta_W,numD] = bu5(noiseRation, varargin )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
%M=15
numD=1;
f=zeros(1,4);
f(1)=0.5;
for i=2:15
    f(i)=0.5-0.025*i;
end
f=f';
bm='independent';
variance=f.*(1-f);
beta_W=1./variance;
end
