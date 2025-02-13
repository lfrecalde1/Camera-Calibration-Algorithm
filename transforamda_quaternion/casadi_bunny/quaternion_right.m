function [output_aux] = quaternion_right(q1)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
qw = q1(1);
qx = q1(2);
qy = q1(3);
qz = q1(4);


Q1 = [qw, -qx, -qy, -qz;...
      qx, qw, qz, -qy;...
      qy, -qz, qw, qx;...
      qz, qy, -qx, qw];
  
output_aux= Q1;
end