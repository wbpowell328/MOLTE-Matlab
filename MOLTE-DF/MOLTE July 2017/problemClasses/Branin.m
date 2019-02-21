function [mu , beta_W, numD]=Branin(noiseRatio, varargin )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
bm='correlated';
numD=2;
xx=-5:1:10;
yy=0:1:15;
[x,y]=meshgrid(xx,yy);
f=fun(x,y);
% surf(x,y,f);
[a,b]=size(f);
mu=max(max(f))-reshape(f,[a*b,1]);
% figure
% plot(mu,'c');
beta_W=1./(noiseRatio*(max(max(f))-min(min(f)))).^2*ones(length(mu),1);
end

function f=fun(x,y)
s=10;
r=6;
a=rand()*2+0.5;
f = a*(y-(5.1/(4*pi^2))*x.^2+5*x/pi-r).^2+s*(1-1/(8*pi))*cos(x)+s;
end
