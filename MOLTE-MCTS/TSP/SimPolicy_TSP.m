function Node = SimPolicy_TSP(rolloutPolicy, Node, root_city, n_cities, costMat)

if length(Node.visited) == n_cities
    Node.V_x = Node.accCost + costMat(Node.city, root_city);
    return;
end

if Node.V_x ~= 0
    return;
end

% Upper bound of costs
ub = sum(sum(costMat));

% Cities that we are solving the TSP for (exclude all previously visited)
cities = [setdiff(1:n_cities, Node.visited), Node.city];

% Current city
curr_city = Node.city;

% Get estimate of best value given the traversed paths
[~, cost] = rolloutPolicy(cities, curr_city, root_city, costMat, ub);

% Add costs up and update value of a node
Node.V_x = Node.accCost + cost;
% fprintf('cost: %d, accumulated cost: %d\n', cost, Node.accCost);
    
end