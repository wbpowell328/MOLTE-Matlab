function [mu , beta_W, numD]=CamelBack(noiseRatio, varargin  )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
bm='correlated';
numD=2;
xx=-2:0.4:2;
yy=-1:0.2:1;
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
f = (4-2.1*x.^2+x.^4/3).*x.^2+x.*y+(-4+4*y.^2).*y.^2+rand()*5;
end
