function f = fitness(index)

global deltas pos_actual pos_lider_actual

aux = squeeze(pos_actual);

if index > 1
    f =  deltas(:, index) - aux(:, index);
else
    f = -pos_lider_actual;
end

end

