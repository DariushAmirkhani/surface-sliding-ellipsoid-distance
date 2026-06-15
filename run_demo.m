% RUN_DEMO
% Simple launcher for the surface-sliding ellipsoid-distance demonstration.
%
% Instructions:
%   1. Open MATLAB.
%   2. Set the current folder to this repository folder.
%   3. Run:
%          run_demo
%
% This will execute the main script:
%   surface_sliding_ellipsoid_distance.m

clear; clc; close all;

repo_dir = fileparts(mfilename('fullpath'));
addpath(repo_dir);

fprintf('Running surface-sliding ellipsoid-distance demo...\n\n');

surface_sliding_ellipsoid_distance;

fprintf('\nDemo finished. Check the MATLAB figures and SurfaceSliding_Log.txt.\n');
