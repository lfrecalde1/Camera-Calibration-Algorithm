clc;close all;clear

load("points.mat");

% Dual quaternion optimization
C1 = zeros(4, 4);
C2 = zeros(4, 4);

npts = size(dest,1);

for i = 1:npts
    b =    source(i,:);
    a =   dest(i,:);

    axbx = a(1) * b(1);
    ayby = a(2) * b(2);
    azbz = a(3) * b(3);
    axby = a(1) * b(2);
    aybx = a(2) * b(1);
    axbz = a(1) * b(3);
    azbx = a(3) * b(1);
    aybz = a(2) * b(3);
    azby = a(3) * b(2);

    C1(1, 1) = C1(1, 1) + axbx - azbz - ayby;
    C1(2, 2) = C1(2, 2) + ayby - azbz - axbx;
    C1(3, 3) = C1(3, 3) + azbz - axbx - ayby;
    C1(4, 4) = C1(4, 4) + axbx + ayby + azbz;
    C1(1, 2) = C1(1, 2) + axby + aybx;
    C1(1, 3) = C1(1, 3) + axbz + azbx;
    C1(1, 4) = C1(1, 4) + aybz - azby;
    C1(2, 3) = C1(2, 3) + azby + aybz;
    C1(2, 4) = C1(2, 4) + azbx - axbz;
    C1(3, 4) = C1(3, 4) + axby - aybx;

    C2(2) = C2(2) + a(3) + b(3);
    C2(3) = C2(3) - a(2) + b(2);
    C2(4) = C2(4) + a(1) - b(1);
    C2(7) = C2(7) + a(1) + b(1);
    C2(8) = C2(8) + a(2) - b(2);
    C2(12) = C2(12) + a(3) - b(3);
  
end

C1(2, 1) = C1(1, 2);
C1(3, 1) = C1(1, 3);
C1(4, 1) = C1(1, 4);
C1(3, 2) = C1(2, 3);
C1(4, 2) = C1(2, 4);
C1(4, 3) = C1(3, 4);
C2(5) = -C2(2);
C2(9) = -C2(3);
C2(13) = -C2(4);
C2(10) = -C2(7);
C2(14) = -C2(8);
C2(15) = -C2(12);

C1 = -2.0 * C1;
C2 = 2.0 * C2;

A = (0.25 / double(npts)) * (C2' * C2) - C1;

% EigenSolver
[evectors, evals] = eig(A);

% Encuentra el índice del máximo eigenvalor
[~, i] = max(diag(evals));

% Extrae el eigenvector asociado al máximo eigenvalor
qmat = evectors(:, i)
%qmat = [qmat(4) qmat(3) qmat(2) qmat(1)]';
%qmat = qmat(end:-1:1);
smat = -(0.5 / double(npts)) * C2 * qmat;

% Convierte a cuaterniones
q = quaternion(qmat(4), qmat(3), qmat(2), qmat(1));
s = quaternion(smat(4), smat(3), smat(2), smat(1));

% Calcula el cuaternión dual
t = s * conj(q);


q_aux = compact(q);

% Convierte el cuaternión a una matriz de rotación
R = quat2rotm([qmat(4), qmat(3), qmat(2) ,qmat(1)]); 

% Asigna los elementos de la matriz de rotación a la matriz de transformación
Rt(1:3, 1:3) = R(1:3, 1:3);

%convierto el quaternion t en un vector
t_aux = compact(t);
tv = [-t_aux(2), -t_aux(3),-t_aux(4)];

% Ajusta las últimas tres columnas de la matriz de transformación
Rt(1, 4) = tv(1);
Rt(2, 4) = tv(2);
Rt(3, 4) = tv(3);

%% transformacion



points_out = dest + tv;

points_out = R * points_out';
points_out = points_out';

figure
plot3(source(:,1),source(:,2),source(:,3),'r.'); hold on
plot3(dest(:,1),dest(:,2),dest(:,3),'b.'); hold on
plot3(points_out(:,1),points_out(:,2),points_out(:,3),'go'); 

axis equal

euler = rotm2eul(R)*180/pi;


numeros = [qmat(1), qmat(2), qmat(3) ,qmat(4)];
combinaciones = [];

% for a = numeros
%     for b = numeros(numeros ~= a)
%         for c = numeros(numeros ~= a & numeros ~= b)
%             d = numeros(numeros ~= a & numeros ~= b & numeros ~= c);
%             combinacion = [a, b, c, d];
%             %combinaciones = [combinaciones; combinacion];
%             R_aux = quat2rotm(combinacion); 
%             euler = rotm2eul(R_aux)*180/pi
% 
%         end
%     end
% end


evals
evectors



