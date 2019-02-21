function [obs, newNode] = generateObservations_TSP(oldNode, newNode, costMat, action, type)

cost = costMat(oldNode.city, action);

if strcmp(type, 'D')
    obs = cost;
    newNode.pmf = containers.Map(cost, 1);
else
    % This check isn't absolutely necessary
    if ~strcmp(type, 'S')
        fprintf('Invalid type %s! Using S (stochastic).', type);
    end
    obs = cost + unique(normrnd(0, 0.1, 20, 1));
    newNode.pmf = containers.Map(obs, ones(length(obs),1)/length(obs));
end

end