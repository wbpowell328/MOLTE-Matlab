%{
Generate N sample realizations for each alternative.
%}
function [samples ] = sampleGen( mu, beta_W, N )

 truth=repmat(mu, [1, N]);
 betaW=repmat(beta_W, [1,N]);
 M=length(mu);
 samples=truth+randn(M,N)./sqrt(betaW);

end

