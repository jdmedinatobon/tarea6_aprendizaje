function phi = phi(fi, fj)

global beta

aux = abs(fi-fj);
% 
% if aux <= beta(1)
%     phi = 1;
% else
%     phi = beta(1)/aux;
% 
% end

phi = (aux <= beta) + (aux > beta).*(beta./aux);
phi(isnan(phi)) = 1;

end

