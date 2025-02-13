%% Code to simulate a camera with internal and external parameters
clc; clear; close all;

% Create camera class (same as your original)
% ---------------------------------------------------------------
% Larger radial and tangential distortion
k1 = -0.2;    % bigger negative => stronger barrel distortion
k2 =  0.01; 
k3 =  0;      
p1 =  0.00;   % non-zero tangential distortion
p2 = -0.00;   % negative sign for demonstration

cam = CentralCamera( ...
    'focal', 0.015, ...             % 15 mm focal length
    'pixel', 10e-6, ...             % 10 micrometers per pixel
    'resolution', [1280 1024], ...  % image size: (nu=1280, nv=1024)
    'centre', [640 512], ...        % principal point
    'k', [k1 k2 k3], ...           % radial distortion
    'p', [p1 p2],...
    'distortion', [k1 k2 k3 p1 p2]);                % tangential distortion

K = cam.K;  % camera intrinsics

% Define chessboard in the XY-plane (Z=0)
% ---------------------------------------------------------------
nx = 4; ny = 5; 
dx = 0.1; dy = 0.1;
[x_vals, y_vals] = meshgrid(0:dx:(nx-1)*dx, 0:dy:(ny-1)*dy);
Z = zeros(size(x_vals));
points3D = [x_vals(:), y_vals(:), Z(:)]';
points3D = points3D(:, 2:end);
points3D_H = [points3D; ones(1,size(points3D,2))];

% Storage arrays
samples = 40;
data_uv = ones(3, size(points3D, 2), samples);  % (2, #points, #samples)
data_xy = ones(3, size(points3D, 2), samples); % (2, #points, #samples)

% Save rotations
R_plane = zeros(3, 3, samples);
t_plane = zeros(3, samples);

%% Aux values for normalized
F_identity = [1, 0, 0;...
              0, 1, 0;...
              0, 0, 1];
Identity = [1, 0, 0, 0;...
            0, 1, 0, 0;...
            0, 0, 1, 0];
% Create a figure for live plotting
figure('Name','Image Plane Animation','NumberTitle','off');
axis([0 cam.nu 0 cam.nv]);   % x from [0..1280], y from [0..1024]
axis ij;                     % flip y-axis to match image coords
hold on; grid on;
xlabel('u (pixels)'); ylabel('v (pixels)');
title('Projected Points in the Image Plane');

for k = 1:samples
    
    % ------- 1) Generate random rotation angles (in degrees) -------
    roll_deg  = -40 + 80*rand(1);  % from -5 to +5
    pitch_deg = -40 + 80*rand(1);  % from -5 to +5
    yaw_deg   = -20 + 40*rand(1);  % from -5 to +5
    
    % Convert chosen angles to radians (or fix some to zero as you did)
    roll  = deg2rad(roll_deg);
    pitch = deg2rad(pitch_deg);
    yaw   = deg2rad(yaw_deg);
    
    % ------- 2) Build the rotation matrix -------
    Rx = [1 0 0; 0 cos(roll) -sin(roll); 0 sin(roll) cos(roll)];
    Ry = [cos(pitch) 0 sin(pitch); 0 1 0; -sin(pitch) 0 cos(pitch)];
    Rz = [cos(yaw) -sin(yaw) 0; sin(yaw) cos(yaw) 0; 0 0 1];
    R_plane(:, :, k) = Rz * Ry * Rx;
    
    % ------- 3) Random translation in Z between [1.5, 2.0] -------
    translationZ = 1.5 + (2.0 - 1.5) * rand(1);
    translationy = -0.4 + (0.1 - (-0.4)) * rand(1,1);
    translationx = -0.4 + (0.1 - (-0.4)) * rand(1,1);
    t_plane(:, k) = [translationx; translationy; translationZ];
    
    % ------- 4) Construct full 4x4 transform H_plane -------
    H_plane = [R_plane(:, :, k), t_plane(:, k); 0 0 0 1];
    
    % ------- 5) Transform the chessboard points into world coords -------
    points3D_I_h = H_plane * points3D_H;  
    points3D_I   = points3D_I_h(1:3, :);
    
    % ------- 6) Project onto image plane using the camera intrinsics -------
    %  or directly with 'cam.C'
    values_normalized = F_identity*Identity*points3D_I_h;
    values_normalized = values_normalized(1:2, :)./values_normalized(3, :);
    radius = vecnorm(values_normalized);
    D = 1 + k1*radius.^2 + k2*radius.^4;
    x_warp = values_normalized.*D;
    x_warp = [x_warp; ones(1, length(x_warp))];


    pixels_aux = cam.K * x_warp;         % 3 x N in homogeneous
    pixels_aux = pixels_aux(1:2,:) ./ pixels_aux(3,:);  % 2 x N
    
    % Store the 2D points if you need them later:
    data_uv(1:2, :,k) = pixels_aux;
    data_xy(1:2, :, k) = points3D(1:2, :);
    
    % ----------- Plot the new frame in the same figure -----------
    cla; % clear previous points
    plot(pixels_aux(1,:), pixels_aux(2,:), 'r.', 'MarkerSize',12);
    
    % Optional text
    text(50, 50, sprintf('Iteration %d', k), 'Color','b','FontSize',12);
    
    drawnow;         % update the figure
    pause(0.01);      % short pause to slow down animation
end


%% Section to compue the homography
X = data_xy;
U = data_uv;

%% Set up optimization solver
addpath('/home/fer/casadi-3.6.7-linux64-matlab2018b');
import casadi.*;
H_homograpy = zeros(3, 3, size(X, 3));
for k=1:size(X,3)
    %% Initial homography using cloosed loop form
    [H_init] = homography_analytical(X(:,:, k), U(:,:,k));

    %% Computing values using non-linear least squares
    [H_refined, resnorm] = estimateHomographyNonlinear(X(1:2,:, k), U(1:2,:,k), H_init);

    %% Computing values using casadi
    [H_refined_casadi] = estimateHomographyCasADi(X(1:2,:, k)', U(1:2,:,k)', H_init);

    %% Normalizing both homographies
    H_refined = H_refined/H_refined(3, 3);
    H_refined_casadi = H_refined_casadi/H_refined_casadi(3, 3);
    H_homograpy(:, :, k) = H_refined_casadi;

    %% Computing the real values and estimation
    U_real = U(1:2, :, k);
    U_estimated = H_refined*X(1:3,:, k);
    U_estimated_casadi = H_refined_casadi*X(1:3,:, k);

end
%% Computing intrinsec parameters
[V] = V_big(H_homograpy);

% Perform Singular Value Decomposition
[U_v, S_v, V_v] = svd(V);

% Extract singular values from the diagonal of S into a vector
sing_vals = diag(S_v);

% Find the index of the smallest singular value
[~, idx] = min(sing_vals);

% The solution is the corresponding column of V
b = V_v(:, idx);
B = [b(1), b(2), b(4);...
     b(2), b(3), b(5);...
     b(4), b(5), b(6);];
%% Computing intrinsec parameters
try
    % Try Cholesky of B
    L = chol(B);
    disp('B is positive definite.');
catch
    % If the above fails, try -B
    disp('B is not positive definite. Checking if -B is positive definite...');
    try
        L = chol(-B);
        disp('B is negative definite.');
    catch
        disp('B is neither positive nor negative definite (indefinite).');
    end
end

A_t = L(3, 3)*(inv(L))';
%% The matrix of intrinsec parameters is defined below
A = (A_t');
A_inv = pinv(A);
K = cam.K

%% Empty values for the estimation
R_estimation = zeros(3, 3, samples);
t_estimation = zeros(3, samples);
%% Section to compute the trasnlation and rotation matrix
for k=1:size(X,3)
    h_0 = H_homograpy(:, 1, k);
    h_1 = H_homograpy(:, 2, k);
    h_2 = H_homograpy(:, 3, k);
    %% Computing lambda
    lambda = 1/(norm(A_inv*h_0));

    %% Computing basis of the rotation matrix
    r_0 = lambda*A_inv*h_0;
    r_1 = lambda*A_inv*h_1;
    r_2 = cross(r_0, r_1);
    R = [r_0, r_1, r_2];
    [U_r, S_r, V_r] = svd(R);
    R_estimation(:, :, k) = (U_r*V_r');

    t_estimation(:, k) = lambda*A_inv*h_2;
end

%% COmputing estimation of the distortion in the lens
uv_c = [A(1,3);A(2, 3)];

D = [];
b = [];
for k=1:size(X,3)
    U_real = U(1:2, :, k);

    T_estimated = [R_estimation(:, :, k), t_estimation(:, k); 0 0 0 1];
    U_estimated = A*Identity*T_estimated*[X(1:2,:, k); zeros(1, size(X, 2)); ones(1, size(X, 2))];
    U_estimated = U_estimated(1:2, :)./U_estimated(3, :);
    
    %% Center real values
    U_real_center = U_real - uv_c;

    center_estimated = U_estimated - uv_c;
    radius_estimated = vecnorm(values_normalized);

    aux_radius_square = U_real_center.*radius_estimated.^2;
    aux_radius_square_square = U_real_center.*radius_estimated.^4;

    aux_radius_square_reshape = reshape(aux_radius_square, 2*size(aux_radius_square, 2), 1);
    aux_radius_square_square_reshape = reshape(aux_radius_square_square, 2*size(aux_radius_square_square, 2), 1);

    aux = [aux_radius_square_reshape, aux_radius_square_square_reshape];

    D = [D;aux];

    %% Compute other side of the equation
    error = U_real - U_estimated;
    error_reshape = reshape(error, 2*size(error, 2), 1);
    b = [b;error_reshape];

end
distortion = pinv(D)*b;

[U_d,S_d,V_d] = svd(D, 'econ');


S_plus = pinv(S_d);        % Since S is square in this 3x2 example ('econ'), S_plus = S^(-1).
D_plus = V_d * S_plus * U_d'; 

% Alternatively, we could have done:
% A_plus = pinv(A);

% Step 3: Solve for x
distortion_2 = D_plus * b;
error_estimation_wo_distortion = [];
error_estimation_distortion = [];
for k=1:size(X,3)
    U_real = U(1:2, :, k);

    T_estimated = [R_estimation(:, :, k), t_estimation(:, k); 0 0 0 1];
    U_estimated = A*Identity*T_estimated*[X(1:2,:, k); zeros(1, size(X, 2)); ones(1, size(X, 2))];
    U_estimated = U_estimated(1:2, :)./U_estimated(3, :);

    error_estimation_wo_distortion = [error_estimation_wo_distortion; (U_real' - U_estimated')];

    values_normalized = F_identity*Identity*T_estimated*[X(1:2,:, k); zeros(1, size(X, 2)); ones(1, size(X, 2))];
    values_normalized = values_normalized(1:2, :)./values_normalized(3, :);
    radius = vecnorm(values_normalized);
    D = 1 + distortion_2(1)*radius.^2 + distortion_2(2)*radius.^4;
    x_warp = values_normalized.*D;
    x_warp = [x_warp; ones(1, length(x_warp))];

    U_improved = A * x_warp;        
    U_improved = U_improved(1:2,:) ./ U_improved(3,:)  ;
    error_estimation_distortion = [error_estimation_distortion; (U_real' - U_improved')];
end

%% Check projection error
norm(error_estimation_distortion)
norm(error_estimation_wo_distortion)


x_init = [A(1,1), A(1, 2), A(2,2), A(1,3), A(2, 3), 0, 0];
%% Section to find the best values based on the casadi optimization solver
a_opt_vec = estimateIntrinsectCasADi(X, U, x_init, R_estimation, t_estimation)
A_optimization = [a_opt_vec(1), a_opt_vec(2), a_opt_vec(4);...
    0, a_opt_vec(3), a_opt_vec(5);...
    0, 0, 1];

%% Verify solution with the error
error_estimation_wo_distortion = [];
error_estimation_distortion = [];
for k=1:size(X,3)
    U_real = U(1:2, :, k);

    T_estimated = [R_estimation(:, :, k), t_estimation(:, k); 0 0 0 1];
    U_estimated = A*Identity*T_estimated*[X(1:2,:, k); zeros(1, size(X, 2)); ones(1, size(X, 2))];
    U_estimated = U_estimated(1:2, :)./U_estimated(3, :);

    error_estimation_wo_distortion = [error_estimation_wo_distortion; (U_real' - U_estimated')];

    values_normalized = F_identity*Identity*T_estimated*[X(1:2,:, k); zeros(1, size(X, 2)); ones(1, size(X, 2))];
    values_normalized = values_normalized(1:2, :)./values_normalized(3, :);
    radius = vecnorm(values_normalized);
    D = 1 + a_opt_vec(6)*radius.^2 + a_opt_vec(7)*radius.^4;
    x_warp = values_normalized.*D;
    x_warp = [x_warp; ones(1, length(x_warp))];

    U_improved = A_optimization * x_warp;        
    U_improved = U_improved(1:2,:) ./ U_improved(3,:)  ;
    error_estimation_distortion = [error_estimation_distortion; (U_real' - U_improved')];
end
norm(error_estimation_distortion)
norm(error_estimation_wo_distortion)