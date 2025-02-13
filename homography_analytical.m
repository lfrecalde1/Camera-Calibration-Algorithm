function [H] = homography_analytical(X, U)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% Normalize the data
[H_X] = normalization_values(X);
[H_U] = normalization_values(U);
X = H_X*X;
U = H_U*U;

A = [];
for k=1:size(U, 2)
    % Split values of the model values
    x_0 = X(1, k);
    y_0 = X(2, k);

    % Split values of the sensor values
    u_0 = U(1, k);
    v_0 = U(2, k);

    aux_a = [-x_0, -y_0, -1, 0, 0, 0, u_0*x_0, u_0*y_0, u_0];
    aux_b = [0, 0, 0, -x_0, y_0, -1, v_0*x_0, v_0*y_0, v_0];

    A = [A; aux_a; aux_b];
end

% Perform Singular Value Decomposition
[U, S, V] = svd(A);

% Extract singular values from the diagonal of S into a vector
sing_vals = diag(S);

% Find the index of the smallest singular value
[~, idx] = min(sing_vals);

% The solution is the corresponding column of V
h_nomr = V(:, idx);
H_norm = reshape(h_nomr, 3, 3)';
H = pinv(H_U)*H_norm*H_X;
end