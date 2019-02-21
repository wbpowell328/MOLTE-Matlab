function [tour, val] = GreedyTSP(cities, curr_city, root, dist, ub)
tour = zeros(length(cities)+1,1);
tour(1) = curr_city;
tour(end) = root;
unvisited_cities = setdiff(cities, curr_city);
val = 0;
ind = 2;
while ~isempty(unvisited_cities)
    [v, i] = min(dist(curr_city,unvisited_cities));
    val = val + v;
    curr_city = unvisited_cities(i);
    tour(ind) = curr_city;
    unvisited_cities = setdiff(unvisited_cities, curr_city);
    ind = ind + 1;
end
val = val + dist(tour(end-1), root);
end