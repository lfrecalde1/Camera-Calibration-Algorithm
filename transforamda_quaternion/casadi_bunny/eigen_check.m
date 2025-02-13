% Example matrix
A = [4, -1; 2, 3];

% Eigen decomposition
[V, D] = eig(A);

% Construct matrices Q and Lambda
Q = V;
Lambda = D;

% Verify Eigen Decomposition
reconstructed_A = Q * Lambda * inv(Q);

% Display results
disp('Original Matrix A:');
disp(A);

disp('Eigenvalues:');
disp(diag(D));

disp('Eigenvectors:');
disp(V);

disp('Matrix Q:');
disp(Q);

disp('Matrix Lambda:');
disp(Lambda);

disp('Reconstructed Matrix A:');
disp(reconstructed_A);