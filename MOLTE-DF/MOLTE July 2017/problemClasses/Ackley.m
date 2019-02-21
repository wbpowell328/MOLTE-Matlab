function [mu , beta_W, numD]= Ackley(noiseRatio, varargin)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
bm='correlated';
numD=2;
xx=-3:0.5:3;
yy=-3:0.5:3;
[x,y]=meshgrid(xx,yy);
aa=rand()*20+10;
bb=rand()*0.1+0.1;
cc=(rand()+1.5)*pi;
f=fun(x,y,aa,bb,cc);
% figure
% surf(x,y,f);
[a,b]=size(f);
mu=max(max(f))-reshape(f,[a*b,1]);
% figure
% plot(mu,'c');
beta_W=1./(noiseRatio*(max(max(f))-min(min(f)))).^2*ones(length(mu),1);
end

function f=fun(x,y,aa,bb,cc)
 f=-aa*exp(-bb*sqrt(1/2*(x.^2+y.^2)))-exp(1/2*(cos(cc*x)+cos(cc*y)));
end
