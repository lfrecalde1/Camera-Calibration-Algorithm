function [H_opt, resnorm] = estimateHomographyNonlinear(pts1, pts2, H_init)
%ESTIMATEHOMOGRAPHYNONLINEAR Refines a 3x3 homography using nonlinear optimization.
%
%   [H_opt, resnorm] = estimateHomographyNonlinear(pts1, pts2, H_init)
%
%   INPUT:
%     pts1  - Nx2 matrix of (x, y) coordinates in source image.
%     pts2  - Nx2 matrix of (x, y) coordinates in target image.
%     H_init - 3x3 initial guess for homography (e.g., from DLT).
%
%   OUTPUT:
%     H_opt  - 3x3 refined homography.
%     resnorm - sum of squared residuals from lsqnonlin.

    % Convert initial H to parameter vector
    params0 = H_init(:);

    % Set up optimization options (Levenberg-Marquardt, e.g.)
    options = optimoptions('lsqnonlin', ...
                           'Algorithm','levenberg-marquardt', ...
                           'Display','iter', ...
                           'MaxIterations', 1000);

    % Define the cost function
    costFunc = @(p) homographyResiduals(p, pts1', pts2');

    % Perform nonlinear least squares
    [params_opt, resnorm] = lsqnonlin(costFunc, params0, [], [], options);

    % Reshape solution back into 3x3
    H_opt = reshape(params_opt, [3, 3]);
end

function residuals = homographyResiduals(params, pts1, pts2)
%HOMOGRAPHYRESIDUALS Returns the vector of reprojection errors
%
%   residuals = homographyResiduals(params, pts1, pts2)
%
%   - params is a 9-element vector representing a 3x3 homography.
%   - pts1, pts2 are Nx2 arrays of corresponding points.

    % Reshape the parameter vector into 3x3 matrix
    H = reshape(params, [3, 3]);

    % Number of correspondences
    N = size(pts1, 1);

    % Initialize residual vector (2D error per correspondence)
    residuals = zeros(2*N, 1);

    for i = 1:N
        % Homogeneous coords of source point
        x1_h = [pts1(i, :), 1]';  % 3x1

        % Project using current H
        x2_est_h = H * x1_h;

        % Convert to inhomogeneous
        x2_est = x2_est_h(1:2) / x2_est_h(3);

        % Actual target point
        x2_actual = pts2(i, :)';  % 2x1

        % Compute residual (difference)
        residuals(2*i-1 : 2*i) = x2_actual - x2_est;
    end
end