function R = rotation_zy(alpha, beta)
% ROTATION_ZY  Rotation matrix R = Rz(alpha) * Ry(beta).
%
% This converts LOCAL coordinates (x',y',z') into GLOBAL (x,y,z):
%   [x; y; z] = R * [x'; y'; z']
%
% INPUTS:
%   alpha - rotation angle around GLOBAL Z axis (radians)
%   beta  - rotation angle around GLOBAL Y axis (radians)
%
% OUTPUT:
%   R     - 3x3 rotation matrix = Rz(alpha) * Ry(beta)

% cos/sin of alpha and beta, to build rotation matrices
ca = cos(alpha); sa = sin(alpha);
cb = cos(beta);  sb = sin(beta);

% Rotation about GLOBAL Z axis by angle alpha:
%   [ ca -sa  0 ]
%   [ sa  ca  0 ]
%   [  0   0  1 ]
Rz = [ ca, -sa, 0;
       sa,  ca, 0;
       0 ,   0, 1];

% Rotation about GLOBAL Y axis by angle beta:
%   [  cb  0  sb ]
%   [   0  1   0 ]
%   [ -sb  0  cb ]
Ry = [ cb, 0, sb;
       0 , 1, 0;
      -sb, 0, cb];

% Combined rotation:
% First apply Ry (rotate around Y), then Rz (around Z):
% Total: R = Rz * Ry
R = Rz * Ry;
end
