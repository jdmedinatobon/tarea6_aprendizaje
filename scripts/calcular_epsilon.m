function epsilon = calcular_epsilon(neighbors)

global beta

n_monos = cellfun('length', neighbors);
epsilon = 0.999./(beta.*n_monos);

end

