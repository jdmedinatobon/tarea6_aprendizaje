function f = fitness(index)

global deltas pos_x_actual pos_x_lider_actual

if index > 1
    f =  deltas(index) - pos_x_actual(index);
else
    f = -pos_x_lider_actual;
end

end

