function decisions = generateDecisions_TSP(oldNode, n_cities)
    decisions = setdiff(1:n_cities, oldNode.visited);
end