function [mu , beta_W, numD]=Easom(noiseRatio, varargin )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
bm='correlated';
numD=2;
xx=1:0.4:5;
yy=1:0.4:5;
[x,y]=meshgrid(xx,yy);
f=fun(x,y);
%surf(x,y,f);
[a,b]=size(f);
mu=max(max(f))-reshape(f,[a*b,1]);
% figure
% plot(mu);
beta_W=1./(noiseRatio*(max(max(f))-min(min(f)))).^2*ones(length(mu),1);
end

function f=fun(x,y)
f =-cos(x).*cos(y).*exp(-(x-pi).^2-(y-pi).^2)+rand()*5;
end
