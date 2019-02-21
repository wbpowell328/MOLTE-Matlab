function [tour, minVal, tElapsed] = TSPDP(cities, curr_city, dist, ub)

vals = java.util.HashMap;
import EqualByValueDoubleSet;

n_cities = length(cities);

finalVals = [];
finalNode = [];
tstart = tic;
for s = 2:n_cities
    % All combinations of size s containing current city
    comb = nchoosek(setdiff(cities, curr_city), s-1);
    comb(:,s) = curr_city;
    [r, c] = size(comb);
    
    finalVals = [];
    % Loop through all combinations of size s and compute subproblem
    % solution
    for k = 1:r
        if c == 2
            vals.put(EqualByValueDoubleSet(comb(k,:),comb(k,1)), dist(curr_city, comb(k,1)));
        elseif c > 2
            S = comb(k,:);
            for i = 1:length(S)
                if S(i) ~= curr_city
                    val = ub*ones(1, length(S));
                    for j = 1:length(S)
                        if S(i) ~= S(j) && S(j) ~= curr_city
                            Sp = setdiff(S, S(i));
                            if any(curr_city == Sp)
                                val(j) = vals.get(EqualByValueDoubleSet(Sp, S(j))) ...
                                    + dist(S(j), S(i));
                                if s == n_cities    % Edge case: add cost of getting back to curr_city to complete tour. Comment out if statement if we don't consider the return trip
                                    val(j) = val(j) + dist(S(i), curr_city);
                                end
                                %                                 fprintf('%s ', string(Sp))
                                %                                 fprintf('i: %d, j: %d, dist: %d, val: %d, valf:%d\n', S(i), S(j), dist(S(j), S(i)), vals.get(EqualByValueDoubleSet(Sp, S(j))), val(j));
                            end
                        end
                    end
                    [v_i, ind] = min(val);
                    vals.put(EqualByValueDoubleSet(S, S(i)), v_i);
                    if s == n_cities
                        finalVals(length(finalVals)+1) = v_i;
                        finalNode(length(finalNode)+1) = S(i);
                    end
                end
            end
            
        end
        
    end
end

% Back pass to extract path
[minVal, minInd] = min(finalVals);
currInd = finalNode(minInd);
tour = zeros(1, length(S));
k = 1;
tour(k) = currInd;
Sp = setdiff(S, currInd);

while length(Sp) > 2
    k = k + 1;
    evals = ub*ones(1,length(Sp));
    inds = zeros(1, length(Sp));
    for i = 1:length(Sp)
        val = vals.get(EqualByValueDoubleSet(Sp, Sp(i))) + dist(Sp(i), tour(k-1));
        if ~ isempty(val)
            evals(i) = val;
            inds(i) = Sp(i);
        end
    end
    [~, minInd] = min(evals);
    currInd = inds(minInd);
    tour(k) = currInd;
    Sp = setdiff(Sp, currInd);
end

tour(n_cities-1) = setdiff(Sp, curr_city);
tour(n_cities)   = curr_city;
tour             = [flip(tour), curr_city];
% minVal           = minVal + dist(tour(length(tour)-1), curr_city);
tElapsed = toc(tstart);
end