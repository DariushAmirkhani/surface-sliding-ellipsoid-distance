% ****************************************************************************************************
% *
% * Developed by:
% * Dariush Amirkhani & Prof. Junfeng Zhang
% *
% * Affiliation:
% * School of Engineering and Computer Science
% * Laurentian University, Sudbury, ON, Canada
% *
% * Contact:
% * damirkhani@laurentian.ca | jzhang@laurentian.ca
% *
% * **************************************************************************************************

% Minimum distance between two ellipsoids via surface-sliding method
% with adaptive step size (rollback mechanism).
%
% Algorithm:
%   1. Step size: delta = eta (initial), constant unless rollback occurs.
%   2. Adaptive delta (rollback): if distance increases, rollback to
%      previous position and halve delta (alternating between particles).
%   3. Stopping criterion: gamma -> 0 (angle between outward normal and
%      the line connecting the two surface points).
%
% Inputs per particle:
%   Geometry: a, b, c (semi-axes)
%   Center position: X0 = [X0; Y0; Z0]
%   Orientation: alpha, beta, Gamma (ZYX Euler angles)
%   Initial surface point: theta_0, phi_0 (parametric angles)
%
% Requires:
%   rotation_zyx.m
%   ellipsoid_param.m
%   ellipsoid_point_local.m
%   rotate_system.m
% =========================================================================
clear; clc;

%% === 1) INPUT PARAMETERS ================================================

a0 = 0.2;  b0 = 0.6;  c0 = 0.4;
eps = 0.1;  % scaling factor for particle #2

% --- Algorithm controls ---
eta      = 0.05;
max_iter = 400;

% --- Particle 1 geometry/orientation/position ---
aI = a0;  bI = b0;  cI = c0;
alphaI = 0.0;  betaI = 0.0;  Gamma1 = pi/180*0;
X0_I = [-1.5; 0; 0];

% --- Particle 2 geometry/orientation/position ---
aJ = a0*eps;  bJ = b0*eps;  cJ = c0*eps;
alphaJ = 0;  betaJ = pi/2;  Gamma2 = 0.0;
X0_J = [1.5; 0; 0];

RY = pi/4;  rotate_system;

% --- Initial surface locations (theta, phi) on each particle ---
thetaI = pi + pi/3;  phiI = pi/3;
thetaJ = -pi/4;      phiJ = pi/2;

% --- Stopping criterion: gamma tolerance (radians) ---
tol_gamma = 1e-4;

% --- Pole clamp for phi ---
eps_phi = 1e-6;

% --- Output file ---
out_file = 'SurfaceSliding_Log.txt';

% --- Figure options ---
label_steps     = [5 10 25];
post_steps_plot = 20;

%% === 2) DERIVED QUANTITIES ==============================================

min_len_I = min([aI, bI, cI]);
min_len_J = min([aJ, bJ, cJ]);

delta_I_init = eta;
delta_J_init = eta;
delta_I = delta_I_init;
delta_J = delta_J_init;

fprintf('=== ALGORITHM PARAMETERS ===\n');
fprintf('eta = %.2f\n', eta);
fprintf('delta_I = eta * min(aI,bI,cI) = %.2f * %.2f = %.4f\n', eta, min_len_I, delta_I);
fprintf('delta_J = eta * min(aJ,bJ,cJ) = %.2f * %.2f = %.4f\n', eta, min_len_J, delta_J);
fprintf('tol_gamma = %.2e rad (%.4f deg)\n', tol_gamma, rad2deg(tol_gamma));
fprintf('============================\n\n');

%% === 3) ROTATIONS AND SURFACES ==========================================

R_I = rotation_zyx(alphaI, betaI, Gamma1);
R_J = rotation_zyx(alphaJ, betaJ, Gamma2);

n = 50;
[ThetaGrid, PhiGrid] = meshgrid(linspace(0, 2*pi, n), linspace(0, pi, n));

[XlocI, YlocI, ZlocI] = ellipsoid_param(aI, bI, cI, ThetaGrid, PhiGrid);
ptsI_global = R_I * [XlocI(:)'; YlocI(:)'; ZlocI(:)'] + X0_I;
XgI = reshape(ptsI_global(1,:), size(XlocI));
YgI = reshape(ptsI_global(2,:), size(XlocI));
ZgI = reshape(ptsI_global(3,:), size(XlocI));

[XlocJ, YlocJ, ZlocJ] = ellipsoid_param(aJ, bJ, cJ, ThetaGrid, PhiGrid);
ptsJ_global = R_J * [XlocJ(:)'; YlocJ(:)'; ZlocJ(:)'] + X0_J;
XgJ = reshape(ptsJ_global(1,:), size(XlocJ));
YgJ = reshape(ptsJ_global(2,:), size(XlocJ));
ZgJ = reshape(ptsJ_global(3,:), size(XlocJ));

%% === 4) 3D PLOT =========================================================

fig1 = figure();
subplot(121);
hold on;

surf(XgI, YgI, ZgI, 'FaceColor', [0.07, 0.62, 1], 'EdgeColor', 'none', ...
    'FaceAlpha', 0.4, 'FaceLighting', 'gouraud', 'AmbientStrength', 0.3, ...
    'DiffuseStrength', 0.8, 'SpecularStrength', 0.9, 'SpecularExponent', 25, ...
    'BackFaceLighting', 'reverselit');

surf(XgJ, YgJ, ZgJ, 'FaceColor', [0.7, 0.45, 0.45], 'EdgeColor', 'none', ...
    'FaceAlpha', 0.4, 'FaceLighting', 'gouraud', 'AmbientStrength', 0.3, ...
    'DiffuseStrength', 0.8, 'SpecularStrength', 0.9, 'SpecularExponent', 25, ...
    'BackFaceLighting', 'reverselit');

L = 0.5;
plot3([1 0 0]*L, [0 0 1]*L, [0 0 0]*L, 'k', 'linewidth', 1);
plot3([0 0]*L, [0 0]*L, [0 1]*L, 'k', 'linewidth', 1);
text(L, 0, 0, 'X');
text(0, L, 0, 'Y');
text(0, 0, L, 'Z');

axis equal; grid on; box on;
view(-45, 20);
camlight('headlight');

%% === 5) OPEN LOG FILE ===================================================

fid = fopen(out_file, 'w');
if fid < 0
    error('Could not open output file: %s', out_file);
end

%% === 6) MAIN ITERATION ==================================================

P_path_I = zeros(3, max_iter+1);
P_path_J = zeros(3, max_iter+1);
theta_hist_I = zeros(max_iter+1, 1);
phi_hist_I   = zeros(max_iter+1, 1);
theta_hist_J = zeros(max_iter+1, 1);
phi_hist_J   = zeros(max_iter+1, 1);

dist_prev   = NaN;
P_I_prev    = NaN(3,1);
P_J_prev    = NaN(3,1);
thetaI_prev = thetaI;  phiI_prev = phiI;
thetaJ_prev = thetaJ;  phiJ_prev = phiJ;

k_final        = 0;
gamma_I        = NaN;
gamma_J        = NaN;
dist           = NaN;
rollback_count = 0;

% Determine which particle's step to halve first on rollback
[~, step_half_indicator] = max([delta_I, delta_J]);

for k = 1:max_iter

    % --- Get local geometry ---
    [x_loc_I, n_loc_I, e_theta_loc_I, e_phi_loc_I, ~, ~] = ...
        ellipsoid_point_local(aI, bI, cI, thetaI, phiI);
    [x_loc_J, n_loc_J, e_theta_loc_J, e_phi_loc_J, ~, ~] = ...
        ellipsoid_point_local(aJ, bJ, cJ, thetaJ, phiJ);

    % --- Transform to global ---
    P_I = R_I * x_loc_I + X0_I;
    P_J = R_J * x_loc_J + X0_J;
    n_I_global = R_I * n_loc_I;  n_I_global = n_I_global / norm(n_I_global);
    n_J_global = R_J * n_loc_J;  n_J_global = n_J_global / norm(n_J_global);

    % --- Store history ---
    P_path_I(:, k) = P_I;
    P_path_J(:, k) = P_J;
    theta_hist_I(k) = thetaI;  phi_hist_I(k) = phiI;
    theta_hist_J(k) = thetaJ;  phi_hist_J(k) = phiJ;

    % --- Compute distance ---
    d_vec = P_J - P_I;
    dist  = norm(d_vec);
    lambda_dir_IJ = d_vec / dist;

    % --- Compute gamma: angle between normal and connecting line ---
    cos_gamma_I = dot(lambda_dir_IJ, n_I_global);
    cos_gamma_J = dot(-lambda_dir_IJ, n_J_global);
    cos_gamma_I = max(-1, min(1, cos_gamma_I));
    cos_gamma_J = max(-1, min(1, cos_gamma_J));
    gamma_I = acos(abs(cos_gamma_I));
    gamma_J = acos(abs(cos_gamma_J));

    % --- Compute displacements ---
    if k == 1
        disp1 = 0.0;  disp2 = 0.0;
    else
        disp1 = norm(P_I - P_I_prev);
        disp2 = norm(P_J - P_J_prev);
    end

    % --- Current eta values ---
    eta_I_current = delta_I / min_len_I;
    eta_J_current = delta_J / min_len_J;

    % --- Write to log ---
    fprintf(fid, '%d\t%.16e\t%.16e\t%.16e\t%.16e\t%.16e\t%.16e\t%.16e\t%.16e\t%.16e\t%.16e\t%.16e\t%.16e\t%.16e\n', ...
        k, thetaI, phiI, thetaJ, phiJ, disp1, disp2, dist, delta_I, delta_J, gamma_I, gamma_J, eta_I_current, eta_J_current);

    % --- Print progress ---
    fprintf('Iter %3d: dist=%.6e, gamma_I=%.4f deg, gamma_J=%.4f deg, delta_I=%.4e\n', ...
        k, dist, rad2deg(gamma_I), rad2deg(gamma_J), delta_I);

    % --- Rollback check: if distance increased, halve delta and rollback ---
    if k > 1 && ~isnan(dist_prev) && dist > dist_prev
        rollback_count = rollback_count + 1;
        fprintf('  >> ROLLBACK #%d: dist increased (%.6e > %.6e)\n', rollback_count, dist, dist_prev);
        fprintf('     delta: %.4e -> %.4e (halved)\n', delta_I, delta_I/2);

        % Restore previous position
        thetaI = thetaI_prev;  phiI = phiI_prev;
        thetaJ = thetaJ_prev;  phiJ = phiJ_prev;

        % Halve step size, alternating between particles
        if step_half_indicator == 1
            delta_I = delta_I / 2;
        else
            delta_J = delta_J / 2;
        end
        step_half_indicator = 3 - step_half_indicator;

        % Restore previous state
        P_I  = P_I_prev;
        P_J  = P_J_prev;
        dist = dist_prev;

        continue;
    end

    % --- Save previous state ---
    P_I_prev    = P_I;
    P_J_prev    = P_J;
    dist_prev   = dist;
    thetaI_prev = thetaI;  phiI_prev = phiI;
    thetaJ_prev = thetaJ;  phiJ_prev = phiJ;

    % --- Compute tangent projections ---
    lambda_dir_JI = -lambda_dir_IJ;

    e_theta_I_global = R_I * e_theta_loc_I;
    e_phi_I_global   = R_I * e_phi_loc_I;
    lambda_theta_I = dot(lambda_dir_IJ, e_theta_I_global);
    lambda_phi_I   = dot(lambda_dir_IJ, e_phi_I_global);

    e_theta_J_global = R_J * e_theta_loc_J;
    e_phi_J_global   = R_J * e_phi_loc_J;
    lambda_theta_J = dot(lambda_dir_JI, e_theta_J_global);
    lambda_phi_J   = dot(lambda_dir_JI, e_phi_J_global);

    % --- Normalize tangent direction ---
    Mi = sqrt(lambda_theta_I^2 + lambda_phi_I^2);
    Mj = sqrt(lambda_theta_J^2 + lambda_phi_J^2);

    % --- Compute parameter updates ---
    if Mi > 1e-15
        d_thetaI = delta_I * (lambda_theta_I / Mi);
        d_phiI   = delta_I * (lambda_phi_I / Mi);
    else
        d_thetaI = 0;  d_phiI = 0;
    end

    if Mj > 1e-15
        d_thetaJ = delta_J * (lambda_theta_J / Mj);
        d_phiJ   = delta_J * (lambda_phi_J / Mj);
    else
        d_thetaJ = 0;  d_phiJ = 0;
    end

    % --- Update angles ---
    thetaI = thetaI + d_thetaI;
    phiI   = phiI   + d_phiI;
    thetaJ = thetaJ + d_thetaJ;
    phiJ   = phiJ   + d_phiJ;

    % --- Clamp phi to avoid poles ---
    phiI = min(max(phiI, eps_phi), pi - eps_phi);
    phiJ = min(max(phiJ, eps_phi), pi - eps_phi);
end

if k_final == 0
    k_final = min(k, max_iter);
    fprintf('\n*** MAX ITERATIONS REACHED (%d)\n', max_iter);
end

fclose(fid);

%% === 7) TRUNCATE HISTORIES ==============================================

P_path_I = P_path_I(:, 1:k_final);
P_path_J = P_path_J(:, 1:k_final);
theta_hist_I = theta_hist_I(1:k_final);
phi_hist_I   = phi_hist_I(1:k_final);
theta_hist_J = theta_hist_J(1:k_final);
phi_hist_J   = phi_hist_J(1:k_final);

%% === 8) FINAL RESULTS ===================================================

[x_loc_I_final, n_loc_I_final, ~, ~, ~, ~] = ellipsoid_point_local(aI, bI, cI, thetaI, phiI);
[x_loc_J_final, n_loc_J_final, ~, ~, ~, ~] = ellipsoid_point_local(aJ, bJ, cJ, thetaJ, phiJ);

P_I_final = R_I * x_loc_I_final + X0_I;
P_J_final = R_J * x_loc_J_final + X0_J;
d_final   = norm(P_J_final - P_I_final);

n_I_final = R_I * n_loc_I_final;  n_I_final = n_I_final / norm(n_I_final);
n_J_final = R_J * n_loc_J_final;  n_J_final = n_J_final / norm(n_J_final);

eta_final_I = delta_I / min_len_I;
eta_final_J = delta_J / min_len_J;

fprintf('\n=== FINAL RESULTS ===\n');
fprintf('Output file: %s\n', out_file);
fprintf('Total iterations: %d\n', k_final);
fprintf('Total rollbacks: %d\n', rollback_count);
fprintf('Final distance: %.16e\n', d_final);
fprintf('Final gamma_I: %.6e rad (%.6f deg)\n', gamma_I, rad2deg(gamma_I));
fprintf('Final gamma_J: %.6e rad (%.6f deg)\n', gamma_J, rad2deg(gamma_J));
fprintf('--- Step size ---\n');
fprintf('Initial delta_I: %.4e  |  Initial eta: %.4f\n', delta_I_init, eta);
fprintf('Final delta_I:   %.4e  |  Final eta:   %.6f\n', delta_I, eta_final_I);
fprintf('Initial delta_J: %.4e  |  Initial eta: %.4f\n', delta_J_init, eta);
fprintf('Final delta_J:   %.4e  |  Final eta:   %.6f\n', delta_J, eta_final_J);
if rollback_count > 0
    fprintf('Delta/eta reduced by factor: 2^%d = %d\n', rollback_count, 2^rollback_count);
end

%% === 9) ADD PATHS TO 3D PLOT ============================================

plot3(P_path_I(1,:), P_path_I(2,:), P_path_I(3,:), 'b', 'LineWidth', 1);
plot3(P_path_J(1,:), P_path_J(2,:), P_path_J(3,:), 'r', 'LineWidth', 1);

ms = 6;

% Start points
plot3(P_path_I(1,1), P_path_I(2,1), P_path_I(3,1), 'o', 'MarkerSize', ms, ...
    'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'k', 'LineWidth', 1);
plot3(P_path_J(1,1), P_path_J(2,1), P_path_J(3,1), 'o', 'MarkerSize', ms, ...
    'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'k', 'LineWidth', 1);

% Final points
plot3(P_I_final(1), P_I_final(2), P_I_final(3), 's', 'MarkerSize', ms, ...
    'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'k', 'LineWidth', 1);
plot3(P_J_final(1), P_J_final(2), P_J_final(3), 's', 'MarkerSize', ms, ...
    'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'k', 'LineWidth', 1);

% Minimum distance line
plot3([P_I_final(1) P_J_final(1)], [P_I_final(2) P_J_final(2)], ...
      [P_I_final(3) P_J_final(3)], 'k--', 'LineWidth', 1);

% Intermediate labeled steps
label_steps = label_steps(label_steps >= 1 & label_steps <= k_final);
for kk = label_steps
    plot3(P_path_I(1,kk), P_path_I(2,kk), P_path_I(3,kk), '^', 'MarkerSize', ms, ...
        'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'k', 'LineWidth', 0.5);
    plot3(P_path_J(1,kk), P_path_J(2,kk), P_path_J(3,kk), '^', 'MarkerSize', ms, ...
        'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'k', 'LineWidth', 0.5);
end

rotate3d on;

%% === 10) CONVERGENCE AND THETA-PHI PLOTS ================================

% --- Distance convergence ---
distances = zeros(k_final, 1);
for ii = 1:k_final
    distances(ii) = norm(P_path_J(:,ii) - P_path_I(:,ii));
end

subplot(222); hold on;
plot(distances, 'k-');
dext = 3 - aI - cJ;
h = plot([0 max_iter], [1 1]*dext, 'b--');
set(h, 'HandleVisibility', 'off');

str = ['Particle #2 scaling factors: [', num2str(aJ/a0), ...
    ',  ', num2str(bJ/b0), ',   ', num2str(cJ/c0), ']'];
legend(str);

% --- Theta-phi polar map ---
subplot(224); hold on;

% Constant-phi rings
a = (0:200)/200 * 2*pi;
ca = cos(a);  sa = sin(a);
for kk = 1:6
    plot(kk*ca, kk*sa, 'k-');
end

% Constant-theta rays
for kk = 0:11
    ang = 2*pi/12 * kk;
    plot([0 cos(ang)*6], [0 sin(ang)*6], 'k-');
end

% Convert (theta, phi) history to polar coordinates
npts = length(theta_hist_I);
xI = zeros(npts, 1);  yI = zeros(npts, 1);
xJ = zeros(npts, 1);  yJ = zeros(npts, 1);
for ii = 1:npts
    r = phi_hist_I(ii) / pi * 6;
    xI(ii) = r * cos(theta_hist_I(ii));
    yI(ii) = r * sin(theta_hist_I(ii));

    r = phi_hist_J(ii) / pi * 6;
    xJ(ii) = r * cos(theta_hist_J(ii));
    yJ(ii) = r * sin(theta_hist_J(ii));
end

plot(xI, yI, 'b-', 'LineWidth', 2);
plot(xJ, yJ, 'r-', 'LineWidth', 2);

ms = 10;

% Start points
plot(xI(1), yI(1), 'bo', 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'k', 'MarkerSize', ms);
plot(xJ(1), yJ(1), 'ro', 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'k', 'MarkerSize', ms);

% Final points
plot(xI(end), yI(end), 'bs', 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'k', 'MarkerSize', ms);
plot(xJ(end), yJ(end), 'rs', 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'k', 'MarkerSize', ms);

axis equal;