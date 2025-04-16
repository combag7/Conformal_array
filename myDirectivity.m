function D_db = myDirectivity(F_norm, phi_vector, theta_vector)
% Calculates directivity using radiation pattern data 

phi_unique = unique(phi_vector);
theta_unique = unique(theta_vector);

F_grid = reshape(F_norm, length(theta_unique), length(phi_unique));


I = cumtrapz(deg2rad(theta_unique), cumtrapz(deg2rad(phi_unique), F_grid.^2 .* sind(theta_unique), 2));   
Power = I(end) / (4 * pi);

% Calculate directivity 
directivity = F_norm.^2/Power;
D_db = changem(10 .* log10(directivity), -40, -inf);
D_db(D_db < -40) = -40;
end