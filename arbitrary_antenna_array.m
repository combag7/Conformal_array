clc; close all; clear all;

addpath('Functions');
%% Setup farfield data

% Reading farfield file
farfield_data1 = table2array(readtable('Import data\my_patch.txt'));

Angles(:, 1) = farfield_data1(:, 1); % Phi [0, 360]
Angles(:, 2) = farfield_data1(:, 2); % Theta [0, 180]

E_theta = farfield_data1(:, 3) + 1j * farfield_data1(:, 4);
E_phi = farfield_data1(:, 5) + 1j * farfield_data1(:, 6); 

% Creating view points array
wave_dir = [sind(Angles(:, 2)) .* cosd(Angles(:, 1)), sind(Angles(:, 2)) .* sind(Angles(:, 1)), cosd(Angles(:, 2))];

% Defining Excitation parameters
Frequency = 14.25e9;
k = 2 * pi * Frequency/physconst('LightSpeed');
Wavelength = physconst('LightSpeed')/Frequency;

%% Creating arrays with elements' positions and arbitrary rotation

num_el = 64;    % Number of elements
len_el = 0.011; % Length of one element (consider them square shaped)

% Spacing megacell
Position = zeros(sqrt(num_el), 3, sqrt(num_el)); % Lists in array = rows in megacell, rows in array = columns in megacell, columns in array = 3 coordinates of each cell
Rotation = zeros(sqrt(num_el), sqrt(num_el));

% Creating square megacell
for i = 2:sqrt(num_el)
    Position(i, 1, :) = Position(i - 1, 1, :) + len_el;
    Position(:, 2, i) = Position(:, 2, i - 1) + len_el;
end

Position = permute(Position, [1 3 2]);
Position = reshape(Position, num_el, 3);

% Rotating elements in megacell
for n = 1:2:sqrt(num_el)
    v = 90 * ones(n, 1);
    vu = diag(v, sqrt(num_el) - n);
    Rotation = Rotation + vu;
    vl = diag(v, -(sqrt(num_el) - n));
    Rotation = Rotation + vl;
end
%Rotation = Rotation - diag(90 * ones(sqrt(nl), 1));

% Rotation flatten
Rotation = reshape(Rotation, [], 1);

dot_kr = zeros(1, length(wave_dir));

%% Beam scanning

% For beam scanning create a linear phase shift mask depending on the tilt
% angle. Apply phase-mask to existing ExcitationPhase array

% Define maximum scanning (tilt) angle in degrees
scan_angle = 0; % 

if scan_angle == 0
    ExcitationPhase = Rotation;
else 
% Calculate linear phase shift 
p = sind(scan_angle) * k;     % radians! per 1 meter

% Creating phase shift array
LinearPhaseShift = zeros(sqrt(num_el), sqrt(num_el));

for m = 1:sqrt(num_el)
    PhaseShift = p * (m - 1) * len_el;

    while abs(PhaseShift) > 2 * pi
        if PhaseShift > 0
        PhaseShift = PhaseShift - 2 * pi;
        elseif PhaseShift < 0
                PhaseShift = PhaseShift + 2 * pi;
        end
    end

    LinearPhaseShift(:, m) = rad2deg(PhaseShift); 
end

% PhaseShift flatten
LinearPhaseShift = reshape(LinearPhaseShift, [], 1);

ExcitationPhase = Rotation + LinearPhaseShift;         % (psi_1, psi_2, ...)
end

% Plotting Excitation phase mask
A = reshape(ExcitationPhase, sqrt(num_el), sqrt(num_el));

figure;
imagesc(A); % Plot heatmap
colormap("parula"); 
colorbar;
caxis([min(A(:)), max(A(:))]);
set(gca, 'YDir', 'normal');
title('Excitation Phase of the array');
axis equal tight;

%% Calculating farfield (Arbitrary array case)

Array = Radiation_Pattern(Frequency, E_phi, E_theta, wave_dir, Position, Rotation, ExcitationPhase);

F_norm_array = Array.F_abs/max(Array.F_abs);

% Total farfield
figure
patternCustom(F_norm_array, Angles(:, 2), Angles(:, 1))
title('Total radiation pattern of array')

%% Calculating directivity
Array.Directivity = myDirectivity(F_norm_array, Angles(:, 1), Angles(:, 2));

disp(['Directivity (array): ', num2str(max(Array.Directivity)), ' dB'])

%% Plotting directivity in polar (Array)
% PatternCustom works in mysterious ways so I want to manually select
% Directivity values and plot them in polar. More freedom of choice!

% Selecting values in slice
slice_val = 90;

indices1 = find(Angles(:, 1) == slice_val);
indices2 = find(Angles(:, 1) == slice_val + 180);

theta1 = Angles(indices1, 2);
dir1 = Array.Directivity(indices1);

theta2 = Angles(indices2, 2);
dir2 = Array.Directivity(indices2);

% Convert theta to radians for polarplot
theta1_rad = deg2rad(theta1);
theta2_rad = deg2rad(theta2);

% Create polar plot
figure;
polarplot(theta1_rad, dir1, 'r-', 'LineWidth', 1.5); hold on;
polarplot(-theta2_rad, dir2, 'r-', 'LineWidth', 1.5);

% Save data to file
save('Exporting data\Directivity_values_Marat.mat', "theta1_rad", "theta2_rad", "dir1", "dir2"); 

% Set radius limits (directivity range)
rlim([-40 25]);

% Set theta axis so 0° is at the top
set(gca, 'ThetaZeroLocation', 'top');
set(gca, 'ThetaDir', 'clockwise');  % Adjust if needed

% Move radial axis to a different location (e.g., 135°)
set(gca, 'RAxisLocation', 135);

% Customize radial tick labels
rticks([-40 -30 -20 -10 0 10]);
% rticklabels({'-40 dB', '-30 dB', '-20 dB', '-10 dB', '0 dB'});

% Labels and legend
legend({'Slice at \phi = 90°', 'Slice at \phi = 270°'}, 'Location', 'south');
title('Directivity Slices at \phi = 90° and 270°');
hold off;

%% Calculating PLF (Array)

% Red Hot Chilli Peppers
PLF_r_array = PLF(Array.E_phi, Array.E_theta, "right");

figure 
patternCustom(PLF_r_array.PLF_pat, Angles(:, 2), Angles(:, 1))
title('Co-Pol (RHCP)')

% LHCP
PLF_l_array = PLF(Array.E_phi, Array.E_theta, "left");

figure 
patternCustom(PLF_l_array.PLF_pat, Angles(:, 2), Angles(:, 1))
title('Cross-Pol (LHCP)')

% Polarization quallity (another way to evaluate quality?)
Pol_Qual = 10 * log10(max(PLF_l_array.PLF_pat)/max(PLF_r_array.PLF_pat));
disp(['Polarization quality (array): ', num2str(Pol_Qual), ' dB'])

%% Calculating directivity (Polarization)

PLF_r_dir = changem(PLF_r_array.PLF, 1, 0) .* Array.Directivity;

PLF_r_dir_max = max(PLF_r_dir);
disp(['Directivity of array (RHCP): ', num2str(max(PLF_r_array.PLF .* Array.Directivity)), ' dB'])

PLF_l_dir = changem(PLF_l_array.PLF, 1, 0) .* Array.Directivity;

disp(['Directivity of array (max, LHCP): ', num2str(max(PLF_l_array.PLF .* Array.Directivity)), ' dB'])

disp(['Directivity of array (observation point, LHCP): ', num2str(PLF_l_dir(find(PLF_r_dir==PLF_r_dir_max, 1))), ' dB'])

%% Plotting polar slice 

dir1_r = PLF_r_dir(indices1);
dir2_r = PLF_r_dir(indices2);

dir1_l = PLF_l_dir(indices1);
dir2_l = PLF_l_dir(indices2);

% Create polar plot
figure;
polarplot(theta1_rad, dir1_r, 'b-', 'LineWidth', 1.5); hold on;
polarplot(-theta2_rad, dir2_r, 'b-', 'LineWidth', 1.5); hold on;
polarplot(theta1_rad, dir1_l, 'r--', 'LineWidth', 1.5); hold on;
polarplot(-theta2_rad, dir2_l, 'r--', 'LineWidth', 1.5);

% Set radius limits (directivity range)
rlim([-40 max([dir1_r; dir2_r])]);

% Set theta axis so 0° is at the top
set(gca, 'ThetaZeroLocation', 'top');
set(gca, 'ThetaDir', 'clockwise');  % Adjust if needed

% Move radial axis to a different location (e.g., 135°)
set(gca, 'RAxisLocation', 135);

% Customize radial tick labels
rticks([-40 -30 -20 -10 0 10]);
% rticklabels({'-40 dB', '-30 dB', '-20 dB', '-10 dB', '0 dB'});

% Labels and legend
legend({'RHCP (\phi = 90°)', 'RHCP (\phi = 270°)', 'LHCP (\phi = 90°)', 'LHCP (\phi = 270°)'}, 'Location', 'south');
title('PLF Directivity Slices at \phi = 90° and 270°');
hold off;
