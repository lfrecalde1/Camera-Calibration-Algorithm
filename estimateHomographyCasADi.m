function H_opt = estimateHomographyCasADi(pts1, pts2, H_init)

addpath('/home/fer/casadi-3.6.7-linux64-matlab2018b');
import casadi.*;

    % Convert to double precision if not already
    pts1 = double(pts1);
    pts2 = double(pts2);

    N = size(pts1,1);
    assert(size(pts2,1) == N, 'pts1 and pts2 must have the same number of points.');

    % Decision variable: 9 parameters for the 3x3 H matrix
    H_sym = SX.sym('H_sym', 9, 1);    % a 9x1 symbolic variable
    H_mat = reshape(H_sym, 3, 3);     % reshape into 3x3

    % Build the sum-of-squared-reprojection-errors cost
    cost = SX(0);
    for i = 1:N
        x1 = pts1(i,1);
        y1 = pts1(i,2);
        x2 = pts2(i,1);
        y2 = pts2(i,2);

        denom   = H_mat(3,1)*x1 + H_mat(3,2)*y1 + H_mat(3,3);
        x2_est  = (H_mat(1,1)*x1 + H_mat(1,2)*y1 + H_mat(1,3)) / denom;
        y2_est  = (H_mat(2,1)*x1 + H_mat(2,2)*y1 + H_mat(2,3)) / denom;

        cost = cost + (x2 - x2_est)^2 + (y2 - y2_est)^2;
    end

    % Define the NLP (no constraints g(x) in this basic example)
    nlp = struct('x', H_sym, 'f', cost);

    % IPOPT solver options (you can tune these as needed)
    opts = struct;
    opts.print_time = 0;
    opts.ipopt.print_level = 3;  % 0..12. Higher -> more verbose

    % Create the CasADi solver using 'ipopt'
    solver = nlpsol('solver', 'ipopt', nlp, opts);

    % Handle initial guess
    if nargin < 3 || isempty(H_init)
        % Default: identity
        H_init_vec = eye(3);
    else
        if isequal(size(H_init), [3,3])
            H_init_vec = H_init;
        elseif isequal(size(H_init), [9,1]) || isequal(size(H_init), [1,9])
            H_init_vec = reshape(H_init,3,3);
        else
            error('H_init must be 3x3 or 9x1.');
        end
    end

    x0 = reshape(H_init_vec, 9, 1);

    % Solve
    sol = solver('x0', x0);

    % Extract the solution
    H_opt_vec = full(sol.x);  % Convert CasADi object to a numeric array
    H_opt = reshape(H_opt_vec, 3, 3);

    % (Optional) Fix scale so that H_opt(3,3) = 1:
    % H_opt = H_opt / H_opt(3,3);

end