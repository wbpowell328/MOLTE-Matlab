function [mu, beta_W, numD]=Griewank(noiseRatio, varargin  )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
bm='correlated';
numD=2;
xx=-3:0.5:3;
yy=-3:0.5:3;
[x,y]=meshgrid(xx,yy);
f=fun(x,y);
%surf(x,y,f);
[a,b]=size(f);
mu=max(max(f))-reshape(f,[a*b,1]);
%figure
%plot(mu);
beta_W=1./(noiseRatio*(max(max(f))-min(min(f)))).^2*ones(length(mu),1);
end

function f=fun(x,y)
 a=rand()*20+30;
 f=1/a*(x.^2+y.^2)-cos(x).*cos(y/sqrt(2))+2;

end