function f = generateGaussianField(x, A, sigma,mean)
    % Generate a 1D Gaussian random field with a given covariance kernel
    
    % Number of points
    n = length(x);      


    % Define the covariance kernel
    cov_kernel = @(x1, x2) A*exp(-(x1 - x2).^2 / (2 * sigma^2));

    % Construct the covariance matrix
    C = zeros(n, n);
    for i = 1:n
        for j = 1:n
            C(i, j) = cov_kernel(x(i), x(j));
        end
    end

    % Generate the Gaussian random field
    L = chol(C + 1e-10*eye(n), 'lower');  % Add small value for numerical stability
    z = randn(n, 1);                      % Standard normal random vector
    f = L * z+mean;                            % Gaussian random
