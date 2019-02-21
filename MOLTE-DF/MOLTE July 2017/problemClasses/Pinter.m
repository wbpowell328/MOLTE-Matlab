function [mu, beta_W, numD] = Pinter(noiseRatio, varargin   )
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
 f=x.^2+2*y.^2+20*sin(y.*sin(x)-x+sin(y)).^2+40*sin(x.*sin(y)-y+sin(x)).^2 +log10(1+(y.^2-2*x+3*y-cos(x)+1).^2)+2*log10(1+2*(x.^2-2*y+3*x-cos(y)+1).^2)+1+rand()*5;
end