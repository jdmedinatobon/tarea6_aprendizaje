function epsilon = calcular_epsilon(neighbors)

global beta_x

num_neighbors = length(neighbors);

epsilon = (0.999/beta_x)*ones(1, num_neighbors);

%Vectorizar esto
for n=1:num_neighbors
    
    n_mono = length(neighbors{n});
    epsilon(n) = epsilon(n)/n_mono;

end

end

