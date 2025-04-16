function result = PLF(E_phi_total, E_theta_total, direction)
% Calculates Polarization Loss Factor

% E-field magnitude
E_abs = sqrt(abs(E_phi_total).^2 + abs(E_theta_total).^2);

if direction == "right"
    % RedHotChilliPeppers
    E_phi_recieve = 1/sqrt(2);
    E_theta_recieve = -1j/sqrt(2);
end
if direction == "left"
    % LHCP
    E_phi_recieve = 1/sqrt(2);
    E_theta_recieve = 1j/sqrt(2);
end

PLF = abs((E_phi_total ./ E_abs * E_phi_recieve) + (E_theta_total ./ E_abs * E_theta_recieve)).^2;

PLF(isnan(PLF)) = 0; 

result.PLF = PLF;
result.PLF_pat = PLF .* E_abs ./ max(E_abs);
end