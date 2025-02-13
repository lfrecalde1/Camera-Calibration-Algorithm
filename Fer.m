%% Code to Calculate the Transformation
clc;clear;close all

%% Load Points of the system
load("points.mat");

%% Number of points
npts = size(dest,1);

%% Data
p = source';
q = dest';

%% Get mean
p_m = (sum(p(1:3,:)')/npts)';
q_m = (sum(q(1:3,:)')/npts)';


%% New Data
pi = p - p_m;
qi = q - q_m;

%% Initial value M
M = zeros(4, 4);
%% Get Sum
for k = 1:length(dest)
   Pi = [0,-pi(1, k), -pi(2, k), -pi(3, k);...
         pi(1, k), 0, pi(3, k),-pi(2, k);...
         pi(2, k),-pi(3, k),0, pi(1, k);...
         pi(3, k),pi(2, k),-pi(1, k),0];
     
  Qi = [0,-qi(1, k), -qi(2, k), -qi(3, k);...
         qi(1, k), 0, -qi(3, k),qi(2, k);...
         qi(2, k),qi(3, k),0,-qi(1, k);...
         qi(3, k),-qi(2, k),qi(1, k),0];
  M = M + Pi'*Qi;   
end
M
% Get Eigen values and vectors
%% Load Casadi to the path of matlab 
addpath('/home/epvs/MATLAB Add-Ons/casadi-3.6.3-linux64-matlab2018b');
import casadi.*;

quaternion = MX.sym('quaternion', 4, 1);
obj = -quaternion'*M*quaternion;
g = [norm(quaternion)]; 
%% General Vector Optimziation Variables
OPT_variables = [quaternion];

nlprob = struct('f', obj, 'x', OPT_variables, 'g', g);
opts = struct;
opts.ipopt.max_iter = 5000;
opts.ipopt.acceptable_tol =1e-20;
opts.ipopt.acceptable_obj_change_tol = 1e-20;

%% Solver of the problem
solver = nlpsol('solver', 'ipopt', nlprob);
args = struct;

args.lbg(1) = 1.0;
args.ubg(1) = 1.0;
% 

%% Initial Conditions System
q_0 = [1; 0; 0; 0];


args.x0 = [q_0]; % initial value of the optimization variables

sol = solver('x0', args.x0, 'lbg', args.lbg(1), 'ubg', args.ubg(1));
Solution = full(sol.x)


quat_opti = Solution
norm(quat_opti)


% Convierte el cuaternión a una matriz de rotación
R = quat2rotm([Solution(1), Solution(2), Solution(3) ,Solution(4)]); 

% Asigna los elementos de la matriz de rotación a la matriz de transformación
Rt(1:3, 1:3) = R(1:3, 1:3);

% %convierto el quaternion t en un vector
% t_aux = compact(t);
% tv = [-t_aux(2), -t_aux(3),-t_aux(4)];

% % Ajusta las últimas tres columnas de la matriz de transformación
% Rt(1, 4) = tv(1);
% Rt(2, 4) = tv(2);
% Rt(3, 4) = tv(3);

 %% transformacion
% 
% points_out = dest + tv;

points_out = R * dest';
points_out = points_out';

figure
plot3(source(:,1),source(:,2),source(:,3),'r.'); hold on
plot3(dest(:,1),dest(:,2),dest(:,3),'b.'); hold on
plot3(points_out(:,1),points_out(:,2),points_out(:,3),'go'); 

axis equal