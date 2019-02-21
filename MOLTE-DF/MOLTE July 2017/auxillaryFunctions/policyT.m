function [mu_est,recommendedArm,  ereward_sum ]=policyT(policy, mu_0,beta_W,covM,samples, mu, alpha)
 [mu_est, choices, recommendedArm]= eval('caller', [policy,'(mu_0,beta_W,covM,samples, alpha, 1)']);
 ereward_sum=sum(choices.*mu);
end
          