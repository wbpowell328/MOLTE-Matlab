function [alpha, newK] = kestens(alphanought,theta, K, prevgradF,gradF)
% Kesten's rule     
% Input 
% alphanought - scaling estimate 
% theta - tunable parameter 
% K - number of times gradient product has changed so far
% prevgradF - gradient from previous iteration
% gradF - gradient from current iteration

% Output 
% alpha - stepsize 
% 
    % compare sign of errors     
%     if (gradF*prevgradF < 0) K = K + 1;
%     end 
    if (dot(gradF, prevgradF) < 0) K = K + 1;
    end 
    
%     alpha = alphanought * (theta)/(theta + K - 1);
    alpha = alphanought * 1 /(theta + K - 1);
    newK = K;
end
