function result = Radiation_Pattern(Frequency, E_phi_element, E_theta_element, wave_dir, Position, Rotation, Excitation_phase)
%RADIATION_PATTERN Calculates radiation pattern of an antenna array

%Radiators are arbitrarily positioned and turned, have excitation phase

% Defining Excitation parameters
k = 2 * pi * Frequency/physconst('LightSpeed');

% Reshaping arrays for circular shift
E_theta_m = reshape(E_theta_element, 181, []); % Make sure to reshape arrays according to number of unique theta values
E_phi_m = reshape(E_phi_element, 181, []);     % Theta is constant within one row, phi is constant within one column

E_phi_total = 0;
E_theta_total = 0; 

dot_kr = zeros(1, length(wave_dir));

for n = 1:numel(Rotation)
    % Circular shift of phi and theta arrays 
    E_theta_s = circshift(E_theta_m, Rotation(n)/1, 2); % The divisor is an increment of phi values
    E_phi_s = circshift(E_phi_m, Rotation(n)/1, 2);
    
    % Flattening arrays 
    E_theta_flat = reshape(E_theta_s, [], 1);
    E_phi_flat = reshape(E_phi_s, [], 1);
    
    % Dot product (k, r_n)
    for i=1:length(wave_dir)
        dot_kr(i) = sum(wave_dir(i,:) .* Position(n,:)); 
    end
    
    E_phi_total = E_phi_total + E_phi_flat .* exp(-1j .* k .* dot_kr.') .* exp(- 1j * deg2rad(Excitation_phase(n))); 
    E_theta_total = E_theta_total + E_theta_flat .* exp(-1j .* k .* dot_kr.') .* exp(- 1j * deg2rad(Excitation_phase(n)));  
end

% Calculating total farfield
F_abs_total = sqrt(abs(E_phi_total).^2 + abs(E_theta_total).^2);

result.F_abs = F_abs_total;
result.E_phi = E_phi_total;
result.E_theta = E_theta_total;
end

