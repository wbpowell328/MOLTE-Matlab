function [f, beta_W,numD] = bu6( noiseRation,varargin )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
%M=20
numD=1;
f=[0.5,0.48,0.37*ones(1,18)]';
bm='independent';
variance=f.*(1-f);
beta_W=1./variance;
end
