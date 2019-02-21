%{
The auxiliary function of AUFs that gives the value of .
%}
function f= AUF_U(p, c, x, meanD, stdD)
%INPUT:
%   x: a vector 
%   p,c parameters of the model
%   meanD:the mean of the Gaussian distribution of D
%   stdD: the standard deviation of D
D=randn(1,1)*stdD+meanD;
f=p*min(x,D)-c*x;
end
