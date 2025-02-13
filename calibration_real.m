%% Code to simulate a camera with internal and external parameters
clc; clear; close all;

load("calibration_values.mat");

%% Section to compue the homography
samples = size(data_xy, 3);
X = data_xy;
U = data_uv;

F_identity = [1, 0, 0;...
              0, 1, 0;...
              0, 0, 1];
Identity = [1, 0, 0, 0;...
            0, 1, 0, 0;...
            0, 0, 1, 0];

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
    radius_estimated = vecnorm(center_estimated);

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


x_init = [A(1,1), A(1, 2), A(2,2), A(1,3), A(2, 3), distortion(1), distortion(2)];
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

%% Map rotations to quaternios
quaternion_estimated = rotm2quat(R_estimation);
X_init = [];
X_init = [X_init;A(1,1); A(1, 2); A(2,2); A(1,3); A(2, 3); ; distortion(1); distortion(2)];
for k=1:size(X, 3)
    x_quaternion = quaternion_estimated(k, 2:4)/quaternion_estimated(k, 1);
    X_init = [X_init;x_quaternion';t_estimation(:, k)];
end

X_opt_vec = cameraCalibrationCasADi(X, U, X_init);

%% Get values from the optimizer
A_final = [X_opt_vec(1), X_opt_vec(2), X_opt_vec(4);...
    0, X_opt_vec(3), X_opt_vec(5);...
    0, 0, 1]
d_final = [X_opt_vec(6), X_opt_vec(7)]
x_quaternions_trans = reshape(X_opt_vec(8:end), 6, samples);
x = x_quaternions_trans(1:3, :);
trans = x_quaternions_trans(4:6, :);
error_estimation_distortion = [];
R_final = zeros(3, 3, samples);
for k=1:size(X, 3)
    U_real = U(1:2, :, k);

    quaternion = [(1 - x(1:3, k)'*x(1:3, k))/(1  + x(1:3, k)'*x(1:3, k)); 2*x(1, k)/(1 + x(1:3, k)'*x(1:3, k)); 2*x(2, k)/(1 + x(1:3, k)'*x(1:3, k)); 2*x(3, k)/(1 + x(1:3, k)'*x(1:3, k))];
    trans_aux = [trans(1, k); trans(2, k); trans(3, k)];
    T_estimated = [quatTorot(quaternion), trans_aux; 0 0 0 1];
    R_final(:, :, k) = quatTorot(quaternion);

    values_normalized = F_identity*Identity*T_estimated*[X(1:2,:, k); zeros(1, size(X, 2)); ones(1, size(X, 2))];
    values_normalized = values_normalized(1:2, :)./values_normalized(3, :);
    radius = vecnorm(values_normalized);
    D = 1 + d_final(1)*radius.^2 + d_final(2)*radius.^4;
    x_warp = values_normalized.*D;
    x_warp = [x_warp; ones(1, length(x_warp))];

    U_improved = A_final * x_warp;        
    U_improved = U_improved(1:2,:) ./ U_improved(3,:)  ;
    error_estimation_distortion = [error_estimation_distortion; norm(U_real' - U_improved')];

end

error_estimation_distortion

figure('Name','Camera Frames in Inertial Frame', 'Position', [100, 100, 1200, 1200]);
hold on; grid on; axis equal;
xlabel('X'); ylabel('Y'); zlabel('Z');
title('Camera Frames in the Inertial Coordinate System');

% Optionally plot the inertial frame axes at the origin:
plot3([0 0.05],[0 0],[0 0],'r','LineWidth',2); % X-axis (red)
plot3([0 0],[0 0.05],[0 0],'g','LineWidth',2); % Y-axis (green)
plot3([0 0],[0 0],[0 0.05],'b','LineWidth',2); % Z-axis (blue)
text(0, 0, 0, 'Inertial', 'FontSize', 8, 'Color', 'k', 'HorizontalAlignment','left');
% Set axis limits: x from -0.3 to 0.3, y from -0.3 to 0.3, and z from 0 to 0.8
xlim([-0.3, 0.3]);
ylim([-0.3, 0.3]);
zlim([0, 0.8]);
% Adjust the view
view(13,20);
% Decide how "long" you want each axis to appear:
scale = 0.03;  % Tweak to suit

% Number of frames:
numFrames = size(R_final,3);

% Create a colormap with a unique color for each frame.
colors = lines(numFrames);  % Alternatively, try jet(numFrames) or hsv(numFrames)

for k = 1:numFrames
    % Translation of frame k (relative to inertial)
    tx = trans(1,k);
    ty = trans(2,k);
    tz = trans(3,k);

    % Rotation matrix for frame k
    Rk = R_final(:,:,k);

    % Extract (scaled) basis vectors for plotting the local axes
    xAxis = Rk(:,1) * scale;
    yAxis = Rk(:,2) * scale;
    zAxis = Rk(:,3) * scale;

    % Build transformation matrix from camera frame to inertial frame
    T_matrix = [Rk, [tx; ty; tz]; 0 0 0 1];

    % Transform points from the camera frame to the inertial frame.
    % Here, X(1:2, :, k) are the 2D coordinates in the camera frame (assumed to lie on Z=0).
    points2D = X(1:2, :, k);
    numPoints = size(points2D, 2);
    % Augment with zeros for Z (since points are on the image plane) and ones for homogeneous coordinates
    points_homogeneous = [points2D; zeros(1, numPoints); ones(1, numPoints)];
    points3D_I = T_matrix * points_homogeneous;  
    
    % Choose a unique color for this frame
    currentColor = colors(k, :);
    % Plot the transformed points in the inertial frame
    plot3(points3D_I(1,:), points3D_I(2,:), points3D_I(3,:), 'o', ...
          'MarkerSize', 2, 'MarkerFaceColor', currentColor);

    % Plot each axis for the camera frame
    % X-axis in red
    plot3([tx, tx + xAxis(1)], [ty, ty + xAxis(2)], [tz, tz + xAxis(3)], 'r', 'LineWidth', 2);
    % Y-axis in green
    plot3([tx, tx + yAxis(1)], [ty, ty + yAxis(2)], [tz, tz + yAxis(3)], 'g', 'LineWidth', 2);
    % Z-axis in blue
    plot3([tx, tx + zAxis(1)], [ty, ty + zAxis(2)], [tz, tz + zAxis(3)], 'b', 'LineWidth', 2);

    % Optionally label the frame
    text(tx, ty, tz, sprintf('Frame %d', k), 'FontSize', 8, 'Color', 'k', 'HorizontalAlignment','left');
    drawnow;         % update the figure
    pause(0.5);      % short pause to slow down animation

end


