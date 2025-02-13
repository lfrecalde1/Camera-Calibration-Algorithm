function [V] = v_qp(H, p, q)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
H_p = H(:, p);
H_q = H(:, q);

H_0p = H_p(1);
H_1p = H_p(2);
H_2p = H_p(3);

H_0q = H_q(1);
H_1q = H_q(2);
H_2q = H_q(3);


V = [H_0p*H_0q;...
     H_0p*H_1q + H_1p*H_0q;...
     H_1p*H_1q;...
     H_2p*H_0q + H_0p*H_2q;...
     H_2p*H_1q + H_1p*H_2q;...
     H_2p*H_2q];
V= V';
end