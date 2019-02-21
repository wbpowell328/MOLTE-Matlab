function [f, beta_W, numD]= bu2(noiseRation, varargin )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
%M=20
numD=1;
f=[0.5,0.42*ones(1,5),0.38*ones(1,14)]';
bm='independent';
variance=f.*(1-f);
beta_W=1./variance;

end
