function phi = phi(fi, fj)

global beta_x 

aux = abs(fi-fj);

if aux <= beta_x
    phi = 1;
else
    phi = beta_x/aux;

end

end

