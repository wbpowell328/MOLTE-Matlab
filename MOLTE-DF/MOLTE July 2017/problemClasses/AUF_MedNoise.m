function  [mu, beta_W, numD]=AUF_MedNoise(varargin)
%In this general AUF setting, we treat the ratio between c and p as an
%unknown parameter coming from a uniform distribution U[0,1]. The
%distribution of D is fixed to be Gaussian distribution with mean 60 and std
%30. The truth is the expected function values over the distribution of the
%ratio which is approximately obtained by sampling num=10000 times.

%if nargin==0, return the expected function value used for prior generation
%and noise level specification
%if nargin==1, return the function value of a sampled parameter form
%U[0,1]as the truth.
numD=1;
noiseRatio=0.3;
p=1;
x=21:120';
num=10000;

meanDr=60;
%stdDr=40;

c=rand(1,1)*p;
meanD=rand()*meanDr;
stdD=noiseRatio*meanD;
cc=zeros(num,length(x)); 

for i=1:num
     cc(i,:) = AUF_U(p,c,x, meanD, stdD)';
end

mu=mean(cc)';
beta_W=1./var(cc)';

end



