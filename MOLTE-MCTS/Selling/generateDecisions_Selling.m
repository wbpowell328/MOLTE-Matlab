function decisions = generateDecisions_Selling(d, inv, t)
leftover = 6*t - inv;
if leftover <= 0
    decisions = 6;
elseif leftover > 6
    decisions = d;
else
    decisions = leftover:6:6;
end
end