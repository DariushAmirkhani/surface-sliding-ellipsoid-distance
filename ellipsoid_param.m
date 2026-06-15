function [x, y, z] = ellipsoid_param(a, b, c, Theta, Phi)
% ELLIPSOID_PARAM  Parametric ellipsoid surface in LOCAL coordinates.
%
% INPUTS:
%   a, b, c  - semi-axes of ellipsoid in local x', y', z' directions
%   Theta    - azimuth angle(s), measured around local z' (0..2π)
%   Phi      - polar angle(s), measured from +z' down to -z' (0..π)
%              (Theta and Phi can be scalars, vectors, or matrices)
%
% OUTPUTS:
%   x, y, z  - coordinates of points on the ellipsoid in LOCAL frame,
%              same size as Theta and Phi.

% x(Theta, Phi) = a * sin(Phi) * cos(Theta)
% y(Theta, Phi) = b * sin(Phi) * sin(Theta)
% z(Theta, Phi) = c * cos(Phi)
% This is the standard parametric formula for an ellipsoid.

x = a .* sin(Phi) .* cos(Theta);   % local x' coordinate
y = b .* sin(Phi) .* sin(Theta);   % local y' coordinate
z = c .* cos(Phi);                 % local z' coordinate
end
