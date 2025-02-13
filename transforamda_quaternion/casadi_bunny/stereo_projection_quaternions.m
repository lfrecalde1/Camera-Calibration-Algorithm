%% Code to Calculate the Transformation
clc;clear;close all

%% Load Points of the system
load("Data.mat");

%% Casadi 
addpath('/home/fer/casadi-3.6.7-linux64-matlab2018b');
import casadi.*;

x = MX.sym('x', 3, 1);
trans = MX.sym('trans', 3, 1);
obj_vect = [];
%% Number of points
npts = size(dest_q,2);

%% Data
p = source_q;
q = dest_q;


%% Get Sum
for k = 1:length(dest_q)
  quaternion = [(1 - x'*x)/(1  + x'*x); 2*x(1)/(1 + x'*x); 2*x(2)/(1 + x'*x); 2*x(3)/(1 + x'*x)];
  
  trans_aux = [trans(1); trans(2); trans(3)];
  costo = quatTorot(quaternion)*p(2:4, k)  + trans_aux - q(2:4, k);
  %costo_n = costo(2:4)'*costo(2:4);
  costo_n = costo(1:3);
  obj_vect = [obj_vect; costo_n];
  
end
obj = obj_vect'*obj_vect;

%% Solve Optimization Problem
% se crea el vector de desiscion solo de una columna
OPT_variables = [x; trans];

nlprob = struct('f', obj, 'x', OPT_variables);

opts = struct;
opts.ipopt.max_iter = 2000;
opts.ipopt.print_level = 5;
opts.print_time = 1;
opts.ipopt.print_frequency_iter = 1;
opts.ipopt.acceptable_tol = 1e-6;
opts.ipopt.acceptable_obj_change_tol = 1e-6;

solver = nlpsol('solver', 'ipopt', nlprob, opts);

args = struct;

args.x0 = [0; 0.5; 0; 0.0; 0.0; 0.0]; % initial value of the optimization variables
tic
sol = solver('x0', args.x0);
solucion = full(sol.x);
toc
tic
sol = solver('x0', args.x0);
solucion = full(sol.x);
toc
r2 = solucion(1:3)'*solucion(1:3);
q0 = (1 - r2)/(1  + r2);
q1 = 2*solucion(1)/(1 + r2);
q2 = 2*solucion(2)/(1 + r2);
q3 = 2*solucion(3)/(1 + r2);

quat_opti = [q0; q1; q2; q3]
trans = [solucion(4:6)]

norm(quat_opti)
points = quatTorot(quat_opti)*source_q(2:4, :) + trans;

plot3(source_q(2,:),source_q(3,:),source_q(4,:),'r.'); hold on
plot3(dest_q(2,:),dest_q(3,:),dest_q(4,:),'b.');
plot3(points(1,:),points(2,:),points(3,:),'mo');
plot3(0,0,0, 'rp')
axis equal;
grid on
grid minor
xlabel('$\textrm{X}[m]$','Interpreter','latex','FontSize',9); ylabel('$\textrm{Y}[m]$','Interpreter','latex','FontSize',9);zlabel('$\textrm{Z}[m]$','Interpreter','latex','FontSize',9);

