function [V] = V_big(H)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
V = [];
for k=1:size(H, 3)
    %% Computing each element of V
    V_12 = v_qp(H(:, :, k), 1, 2);
    V_11 = v_qp(H(:, :, k), 1, 1);
    V_22 = v_qp(H(:, :, k), 2, 2);
    V_aux = [V_12;...
        V_11-V_22];

    %% Stacking the element for each measurement
    V = [V;V_aux];
    
end
end