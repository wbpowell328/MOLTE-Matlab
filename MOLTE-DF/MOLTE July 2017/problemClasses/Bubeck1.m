function [f, beta_W, numD] = bu1(noiseRation,varargin  )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
%M=20 Bernolli
numD=1;
f=[0.5,0.4*ones(1,19)]';
bm='independent';
variance=f.*(1-f);
beta_W=1./variance;
end

