function [theta, choices] = EpsilonGreedyPricing(X,theta_0,MVar,truth,nTrails,C,~,epsilon)
% Weigher: epsilon here
maxTrails = round(epsilon*nTrails);
[m, ~] = size(X);
theta = theta_0;

for i = 1:nTrails
    if i < maxTrails % Random sample
        action = 1+round((m-1)*rand(1));
        choices(i) = action;
    else % Choose max
        [~,action] = max(X*theta);
        choices(i) = action;
    end
    y = truth(action) + sqrt(MVar) *rand(1); % Sample observation
    [B, theta] = update(action, theta, y, C/MVar); % Update matrices
    C = B*MVar;
end

    function [B, theta] = update(action, theta, w, B)
        x = [1; action; action^2];
        e       = w - theta'*x;
        gamma   = 1 + x'*B*x;
        
        theta   = theta + 1/gamma*B*x*e;
        B       = B - 1/gamma*B*(x*x')*B;
    end
end