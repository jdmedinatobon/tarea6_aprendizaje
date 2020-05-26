function f = fitness(index)

global deltas pos_x_actual

f =  deltas(index) - pos_x_actual(index);

end

