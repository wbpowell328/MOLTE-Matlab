
% generate a prior for GPR. If nargin == 0, use default prior;
% otherwise, use the input to construct a prior
function [mu_0, covM, beta_W]=Prior_GPR(varargin)
noiselevel=0.1;

if nargin ==0 %use the default prior
    sigma=100;
    beta=0.01;
    M=100;
else
    paras=varargin{1};
    before=find(paras==',');
    sigma=str2double(paras(1:before-1));
    beforesim=find(paras==';');
    beta=str2double(paras(before+1:beforesim-1));
    M=str2double(paras(beforesim+1:length(paras)));
end

    beta_W=ones(M,1)./(36*sigma*noiselevel^2);  
    mu_0=randn(M,1)*sqrt(sigma);
    covM=zeros(M,M);
    for i=1:M %construct the prior covariance matrix
        for j=1:M
            covM(i,j)=sigma*exp(-beta*(i-j)^2 );
        end
    end

end

