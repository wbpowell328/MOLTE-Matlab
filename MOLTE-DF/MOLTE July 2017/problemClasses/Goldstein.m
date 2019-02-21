function [mu , beta_W, numD]= Goldstein(noiseRatio, varargin )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
bm='correlated';
numD=2;
xx=-3:0.5:3;
yy=-3:0.5:3;
[x,y]=meshgrid(xx,yy);
f=fun(x,y);
% surf(x,y,f);
[a,b]=size(f);
mu=max(max(f))-reshape(f,[a*b,1]);
%figure
%plot(mu);
beta_W=1./(noiseRatio*(max(max(f))-min(min(f)))).^2*ones(length(mu),1);
end

function f=fun(x,y)
 f=(1+(x+y+1).^2.*(19-14*x+3*x.^2-14*y+6*x.*y+3*y.*y)).*(30+(2*x-3*y).^2.*(18-32*x+12*x.^2+48*y-36*x.*y+27*y.^2))++rand()*5;
end

