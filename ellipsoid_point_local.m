function [x_local, n_local, e_theta, e_phi, len_theta, len_phi] = ...
    ellipsoid_point_local(a, b, c, theta, phi)
% ELLIPSOID_POINT_LOCAL
% For a given (theta, phi) on an ellipsoid (in LOCAL frame), compute:
%   - point position x_local
%   - outward normal n_local (unit vector)
%   - tangent direction when theta increases:   e_theta (unit vector)
%   - tangent direction when phi increases:     e_phi   (unit vector)
%   - their lengths before normalization: len_theta = ||∂r/∂theta||,
%                                          len_phi   = ||∂r/∂phi||
%
% INPUTS:
%   a, b, c  - semi-axes of ellipsoid in local x', y', z'
%   theta    - azimuth angle (around local z', 0..2π)
%   phi      - polar angle   (from +z' down to -z', 0..π)
%
% OUTPUTS:
%   x_local   - 3x1 local position vector of the surface point
%   n_local   - 3x1 unit normal at that point (outward)
%   e_theta   - 3x1 unit tangent for increasing theta
%   e_phi     - 3x1 unit tangent for increasing phi
%   len_theta - scalar, length of ∂r/∂theta (before normalization)
%   len_phi   - scalar, length of ∂r/∂phi   (before normalization)

% --------------------------------------------------------------
% 1) POSITION ON THE SURFACE IN LOCAL COORDINATES
% --------------------------------------------------------------
% Same formulas as in ellipsoid_param, but for a single (theta,phi).
x = a * sin(phi) * cos(theta);   % local x' coordinate
y = b * sin(phi) * sin(theta);   % local y' coordinate
z = c * cos(phi);                % local z' coordinate

% Pack into a 3x1 column vector
x_local = [x; y; z];

% --------------------------------------------------------------
% 2) OUTWARD NORMAL VIA GRADIENT OF IMPLICIT FORM
% --------------------------------------------------------------
% Ellipsoid implicit equation:
%   F(x,y,z) = x^2/a^2 + y^2/b^2 + z^2/c^2 - 1 = 0
%
% Gradient ∇F = (2x/a^2, 2y/b^2, 2z/c^2) points outward.
nx = 2*x/(a^2);
ny = 2*y/(b^2);
nz = 2*z/(c^2);
n_vec = [nx; ny; nz];

% Normalize to get a unit normal vector
n_local = n_vec / norm(n_vec);

% --------------------------------------------------------------
% 3) TANGENT VECTORS ∂r/∂theta AND ∂r/∂phi
% --------------------------------------------------------------
% Parametric form r(theta,phi):
%   r_x = a sin(phi) cos(theta)
%   r_y = b sin(phi) sin(theta)
%   r_z = c cos(phi)
%
% Partial derivative wrt theta:
%   ∂r/∂theta = r_theta
r_theta = [-a * sin(phi) * sin(theta);   % d/dtheta of x
            b * sin(phi) * cos(theta);   % d/dtheta of y
            0];                          % z does not depend on theta

% Partial derivative wrt phi:
%   ∂r/∂phi = r_phi
r_phi   = [ a * cos(phi) * cos(theta);   % d/dphi of x
            b * cos(phi) * sin(theta);   % d/dphi of y
           -c * sin(phi)];               % d/dphi of z

% Lengths of these tangent vectors (magnitudes)
len_theta = norm(r_theta);
len_phi   = norm(r_phi);

% --------------------------------------------------------------
% 4) TURN TANGENT VECTORS INTO UNIT VECTORS
% --------------------------------------------------------------
% We want unit directions e_theta, e_phi for the professor's formulas.
% At poles, len_theta or len_phi may be ~0, so we handle that safely.

if len_theta > 1e-12
    e_theta = r_theta / len_theta;   % unit tangent for theta direction
else
    e_theta = [0; 0; 0];             % degenerate case near pole
end

if len_phi > 1e-12
    e_phi = r_phi / len_phi;         % unit tangent for phi direction
else
    e_phi = [0; 0; 0];
end
end
