function [f, beta_W,numD]= bu4(noiseRation,varargin  )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
%M=6
numD=1;
f=[0.5,0.42,0.4,0.4,0.35,0.35]';
bm='independent';
variance=f.*(1-f);
beta_W=1./variance;
end
