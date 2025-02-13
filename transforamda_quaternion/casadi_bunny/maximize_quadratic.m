clc, clear all, close all;
% Example quadratic function coefficients
A = [2, 1; 1, 3];
b = [1; -1];
c = 2;

% Eigen decomposition
[V, D] = eig(A);

% Identify maximum eigenvalue and its eigenvector
[max_eigenvalue, max_eigenvalue_index] = max(diag(D));
max_eigenvector = V(:, max_eigenvalue_index);

% Compute x_max that maximizes the quadratic form
x_max = V * max_eigenvector;

% Evaluate the maximum value of the quadratic function
max_value = quadratic_function(x_max, A, b, c);

% Display results
disp('Original Quadratic Function Coefficients:');
disp('A:');
disp(A);
disp('b:');
disp(b);
disp('c:');
disp(c);

disp('Eigen Decomposition:');
disp('V (Eigenvectors):');
disp(V);
disp('D (Eigenvalues):');
disp(D);

disp('Maximum Eigenvalue and Eigenvector:');
disp('max_eigenvalue:');
disp(max_eigenvalue);
disp('max_eigenvector:');
disp(max_eigenvector);

disp('Point x_max that Maximizes the Quadratic Form:');
disp('x_max:');
disp(x_max);

disp('Maximum Value of the Quadratic Function:');
disp('max_value:');
disp(max_value);

% Quadratic function
function value = quadratic_function(x, A, b, c)
    value = x' * A * x + b' * x + c;
end