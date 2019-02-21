function [mu , beta_W, numD]= APhyper_ellipsoid(noiseRatio, varargin)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
bm='correlated';
numD=2;
xx=-3:0.5:3;
yy=-3:0.5:3;
[x,y]=meshgrid(xx,yy);
% surf(x,y,f);
aa=rand()+1.5
f=fun(x,y,aa);
[a,b]=size(f);
mu=max(max(f))-reshape(f,[a*b,1]);
% figure
% plot(mu,'r');
beta_W=1./(noiseRatio*(max(max(f))-min(min(f)))).^2*ones(length(mu),1);
end

function f=fun(x,y,aa)
 f=x.^2+aa*y.^2;
end
