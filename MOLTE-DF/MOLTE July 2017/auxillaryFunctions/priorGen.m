function  [mu_0, covM, beta_W] = priorGen(problemClass, noiseRatio, Prior, paras)
%generate prior from a given prior or the default value  
if strcmp(Prior, 'Default') % if default
  [mu_0, covM, beta_W]=eval('caller',[['Prior_',problemClass], '()']);
else  %if specify the prior
    if length(paras)~=0  % if specify the parameters in the spreadsheet
         [mu_0, covM, beta_W]=eval('caller',[['Prior_',problemClass], '(paras)']);
    else % use .mat to give a specified prior
        priorMatrix=load(['Prior_', problemClass,noiseRatio,'.mat']);
        mu_0=priorMatrix.mu_0;
        covM=priorMatrix.covM;
        beta_W=priorMatrix.beta_W;
    end
end

end
