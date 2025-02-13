%% Code to Calculate the Transformation
clc;clear;close all

%% Load Points of the system
load("Data.mat");

%% Casadi 
addpath('/home/epvs/MATLAB Add-Ons/casadi-3.6.3-linux64-matlab2018b');
import casadi.*;

quaternion = MX.sym('quaternion', 4, 1);
obj = 0;
%% Number of points
npts = size(dest_q,2);

%% Data
p = source_q;
q = dest_q;

%% Get mean
p_m = (sum(p(1:4,:),2)/npts);
q_m = (sum(q(1:4,:),2)/npts);


%% New Data
pi = p - p_m;
qi = q - q_m;

%% Initial value M
M = zeros(4, 4);
%% Get Sum
for k = 1:npts
   Pi = [pi(1,k),-pi(2,k), -pi(3,k), -pi(4,k);...
         pi(2,k), pi(1,k), pi(4,k),-pi(3,k);...
         pi(3,k),-pi(4,k),pi(1,k), pi(2,k);...
         pi(4,k),pi(3,k),-pi(2,k),pi(1,k)];
     
  Qi = [qi(1,k),-qi(2,k), -qi(3,k), -qi(4,k);...
         qi(2,k), qi(1,k), -qi(4,k),qi(3,k);...
         qi(3,k),qi(4,k),qi(1,k),-qi(2,k);...
         qi(4,k),-qi(3,k),qi(2,k),qi(1,k)];
     
  M = M + Pi'*Qi; 

  %Other Method
  obj = obj + (-quaternion'*(Pi'*Qi)*quaternion);
  
end
%% Solution
[V,D] = eig(M)
% Get Eigen values and vectors
%% Get the maximun Eig values and the Eigvector is the quaternion
[~, i] = max(diag(D));
quat_aux = (V(:, i));
quat_aux_c = congujate_quaternion(quat_aux)
%% Tranform the points to the Origin Frame
dest_q_comp = quaternion_left(quat_aux)*quaternion_right(quat_aux_c)*source_q;

%% Translation 

b = quaternion_left(quat_aux)*quaternion_right(quat_aux_c)*p_m-q_m;
dest_q_comp = dest_q_comp - b;

plot3(source_q(2,:),source_q(3,:),source_q(4,:),'r.'); hold on
plot3(dest_q(2,:),dest_q(3,:),dest_q(4,:),'b.');
plot3(dest_q_comp(2,:),dest_q_comp(3,:),dest_q_comp(4,:),'mo');
plot3(0,0,0, 'rp')
grid on
grid minor
axis equal
xlabel('$\textrm{X}[m]$','Interpreter','latex','FontSize',9); ylabel('$\textrm{Y}[m]$','Interpreter','latex','FontSize',9);zlabel('$\textrm{Z}[m]$','Interpreter','latex','FontSize',9);