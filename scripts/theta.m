function theta = theta(i, j, fi, fj)

global beta pos_actual

% if fi > fj
%     theta = min(pos_actual(j), beta(1));
% elseif fi < fj
%     theta = min(pos_actual(i), beta(1));
% else
%     theta = beta(1);
% end

theta = (fi > fj).*min(pos_actual(:,j), beta) + ...
        (fi < fj).*min(pos_actual(:,i), beta) + ...
        (fi == fj).*beta;

end

