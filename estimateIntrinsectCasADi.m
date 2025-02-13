function a_opt_vec = estimateIntrinsectCasADi(pts1, pts2, a_init, R_estimation, t_estimation)

addpath('/home/fer/casadi-3.6.7-linux64-matlab2018b');
import casadi.*;

% Convert to double precision if not already
pts1 = double(pts1);
pts2 = double(pts2);

N = size(pts1,3);

% Decision variable: 9 parameters for the 3x3 H matrix
a_vector = SX.sym('intrinsect', 7, 1);    % a 9x1 symbolic variable

% Build the sum-of-squared-reprojection-errors cost
cost = 0;
%% Matrix of intrinsect parameters
A = [a_vector(1), a_vector(2), a_vector(4);...
    0, a_vector(3), a_vector(5);...
    0, 0, 1];
%% Aux values for normalized
F_identity = [1, 0, 0;...
    0, 1, 0;...
    0, 0, 1];
Identity = [1, 0, 0, 0;...
    0, 1, 0, 0;...
    0, 0, 1, 0];


for k = 1:N
    % REAL VALUES
    U_real = pts2(1:2, :, k);

    % TRANSFORMATION
    T_estimated = [R_estimation(:, :, k), t_estimation(:, k); 0 0 0 1];

    % NORMALIZED VALUES
    values_normalized = F_identity*Identity*T_estimated*[pts1(1:2,:, k); zeros(1, size(pts1, 2)); ones(1, size(pts1, 2))];
    values_normalized_aux = values_normalized(1:2, :)./values_normalized(3, :);

    % APLYING BARREL DISTORTION
    radius = vecnorm(values_normalized_aux);
    D = 1 + a_vector(6)*radius.^2 + a_vector(7)*radius.^4;
    D_aux = [D;D];
    x_warp = values_normalized_aux.*D_aux;
    x_warp_aux = [x_warp; ones(1, length(x_warp))];


    U_improved = A * x_warp_aux;
    U_normalized_aux = [U_improved(3,:); U_improved(3,:)];
    U_improved_final = U_improved(1:2,:) ./ U_normalized_aux;

    error = (U_real - U_improved_final);
    
    error_reshape = reshape(error, 2*size(error,2), 1);

    cost = cost + error_reshape'*error_reshape;
end
% Define the NLP (no constraints g(x) in this basic example)
nlp = struct('x', a_vector, 'f', cost);

% IPOPT solver options (you can tune these as needed)
opts = struct;
opts.print_time = 0;
opts.ipopt.print_level = 3;  % 0..12. Higher -> more verbose

% Create the CasADi solver using 'ipopt'
solver = nlpsol('solver', 'ipopt', nlp, opts);

x0 = a_init;

% Solve
sol = solver('x0', x0);

% Extract the solution
a_opt_vec = full(sol.x);  % Convert CasADi object to a numeric array
end