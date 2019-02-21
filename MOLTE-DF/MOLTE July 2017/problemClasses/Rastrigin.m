function [mu , beta_W, numD]= Rastrigin(noiseRatio, varargin  )
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
% figure
% plot(mu);
beta_W=1./(noiseRatio*(max(max(f))-min(min(f)))).^2*ones(length(mu),1);
end

function f=fun(x,y)
 a=rand()*10+5;
 f=2*a+(x.^2-10*cos(2*pi*x))+(y.^2-a*cos(2*pi*y));
end
