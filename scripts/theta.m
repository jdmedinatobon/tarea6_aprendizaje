function theta = theta(i, j, fi, fj)

global beta_x pos_x_actual

if fi > fj
    theta = min(pos_x_actual(j), beta_x);
elseif fi < fj
    theta = min(pos_x_actual(i), beta_x);
else
    theta = beta_x;

end

end

