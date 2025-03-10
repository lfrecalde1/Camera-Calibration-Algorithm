clc; clear; close all;

%% Load Points of the system
load("Data.mat");

%% Casadi 
addpath('/home/fer/casadi-3.6.7-linux64-matlab2018b');
import casadi.*;

%% Decision variables
quaternion = MX.sym('quaternion', 4, 1);
trans      = MX.sym('trans', 3, 1);

%% Objective
obj_vect = [];

%% Number of points
npts = size(dest_q,2);

%% Data
p = source_q;
q = dest_q;

%% Compute means
p_m = sum(p, 2)/npts;
q_m = sum(q, 2)/npts;

%% Center data (remove means)
pi = p - p_m;
qi = q - q_m;

%% Build objective
for k = 1:npts
    % Conjugate of quaternion (user-defined function)
    quaternion_c = congujate_quaternion(quaternion);
    
    % Turn translation into a 'quaternion' style vector
    trans_aux = [0; trans(1); trans(2); trans(3)];
    
    % Residual of transformation
    costo   = quaternion_left(quaternion)*quaternion_right(quaternion_c)*p(:, k) + trans_aux - q(:, k);
    
    % Pick whichever norm you want (L2 or squared L2):
    % Option 1: L2 norm (non-smooth but often OK):
    %   costo_n = norm(costo(2:4));  
    % Option 2: squared L2 norm (smooth):
    %   costo_n = dot(costo(2:4), costo(2:4));
    %
    % Example with L2 norm:
    costo_n = costo(2:4);
  obj_vect = [obj_vect; costo_n];
end
obj = obj_vect'*obj_vect;
%% Norm constraint on the quaternion
% Replace "g = 1 - norm(quaternion);" with the algebraic constraint below:
g = norm(quaternion) - 1;

%% Formulate the NLP
OPT_variables = [quaternion; trans];

nlprob = struct('f', obj, 'x', OPT_variables, 'g', g);

%% Solver options
opts = struct;
opts.ipopt.max_iter = 2000;
opts.ipopt.print_level = 5;
opts.print_time = 1;
opts.ipopt.print_frequency_iter = 1;
opts.ipopt.acceptable_tol = 1e-6;
opts.ipopt.acceptable_obj_change_tol = 1e-6;

solver = nlpsol('solver', 'ipopt', nlprob, opts);

%% Constraint bounds (enforce g = 0)
args.lbg(1) = 0.00;
args.ubg(1) = 0.00;

%% Initial guess
args.x0 = [0.00; 0.5; 0; 0; 0.0; 0.0; 0.0];
tic
%% Solve
sol = solver('x0', args.x0,'lbg', args.lbg, 'ubg',args.ubg);
solucion  = full(sol.x);
toc
%% Extract solution
quat_opti = solucion(1:4)
trans_opti = [0; solucion(5:7)]

%% Check the quaternion norm
disp("Optimized quaternion:")
disp(quat_opti)
disp("Norm of quaternion:")
disp(norm(quat_opti))

%% Apply transformation to source data
quat_opti_c = congujate_quaternion(quat_opti);
points = quaternion_left(quat_opti)*quaternion_right(quat_opti_c)*source_q + trans_opti;

%% Plot
figure;
plot3(source_q(2,:), source_q(3,:), source_q(4,:), 'r.'); hold on
plot3(dest_q(2,:),   dest_q(3,:),   dest_q(4,:),   'b.');
plot3(points(2,:),   points(3,:),   points(4,:),   'mo');
plot3(0, 0, 0, 'rp')
axis equal;
grid on
grid minor
xlabel('$\textrm{X}[m]$','Interpreter','latex','FontSize',9);
ylabel('$\textrm{Y}[m]$','Interpreter','latex','FontSize',9);
zlabel('$\textrm{Z}[m]$','Interpreter','latex','FontSize',9);