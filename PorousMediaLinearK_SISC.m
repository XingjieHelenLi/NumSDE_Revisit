% function PorousMediaLinearK_SISC
% Parameters
n = 1001; %2001;                     % Number of grid points
L = 200;                       % Length of the domain [0, L]
x = linspace(0, L, n);        % 1D mesh grid
h = L/(n-1);                  % Grid spacing

% Permeability k(x) as a linear function + random Bernoulli perturbation
k_0 = 1;                      % Permeability constant term
k_1 = 9;                   % Permeability linear term
c = 0; 0.2;                      % Small constant for random perturbation

% Random Bernoulli variable (+1 or -1) for each grid point
random_perturbation = (2 * randi([0, 1], 1, n) - 1);  % Generates +1 or -1

% Permeability function with random perturbations
k = k_0 + k_1 * x + c * random_perturbation;  % k(x) = linear + perturbation

% Define the sigma variable
b = 10.0; 
power = 1.0;
D0 = 1; 
omega = (k/b).^power;
diffusion = sqrt(D0*omega); 

% Source term f(x) as a Gaussian bump
A = -0.01;                        % Amplitude of the Gaussian bump
x_0 = L/2;                    % Center of the Gaussian bump
sigma =L/10;                % Standard deviation of the Gaussian bump
f = A * exp(-(x - x_0).^2 / (2 * sigma^2));  % Gaussian source

% Setup finite difference matrix
K = zeros(n,n);               % Stiffness matrix
F = f(:);                     % Load vector (source term f)

% Loop over grid points to assemble the stiffness matrix
for i = 2:n-1
    % Local diffusion coefficient k at midpoints
    k_left = (k(i) + k(i-1)) / 2;
    k_right = (k(i) + k(i+1)) / 2;
    
    % Finite difference approximation of div(k grad p)
    K(i,i-1) = -k_left / h^2;
    K(i,i) = (k_left + k_right) / h^2;
    K(i,i+1) = -k_right / h^2;
end

% Apply Dirichlet boundary conditions (p(0) = 0 and p(L) = 0)
K(1,:) = 0; K(1,1) = 1;  % p(0) = 0
K(n,:) = 0; K(n,n) = 1;  % p(L) = 0
F(1) = 0; F(n) = 0;

% Solve the linear system
p = K \ F;

% Compute the gradient of p (finite differences)
grad_p = zeros(1,n);  % Gradient of p
for i = 2:n-1
    grad_p(i) = (p(i+1) - p(i-1)) / (2*h);  % Central difference
end
grad_p(1) = (p(2) - p(1)) / h;  % Forward difference at the left boundary
grad_p(n) = (p(n) - p(n-1)) / h;  % Backward difference at the right boundary

% Compute the flux: q(x) = -k(x) * grad p(x)
flux = -k .* grad_p;

% Plot the solution p(x)
figure;

subplot(4,1,1);
plot(x, p, 'b', 'LineWidth', 1.5);
xlabel('x');
ylabel('p(x)');
title('Solution to \nabla \cdot (k(x) \nabla p(x)) = f(x)');
grid on;

% Plot the Gaussian bump (source term)
subplot(4,1,2);
plot(x, f, 'r', 'LineWidth', 1.5);
xlabel('x');
ylabel('f(x)');
title('Gaussian Source Term f(x)');
grid on;

% Plot the flux q(x) = -k(x) * grad p(x)
subplot(4,1,3);
plot(x, flux, 'g', 'LineWidth', 1.5);
xlabel('x');
ylabel('q(x)');
title('Flux: q(x) = -k(x) \cdot \nabla p(x)');
grid on;

% Plot the diffusion coefficient sigma(x)
subplot(4,1,4);
plot(x, diffusion, 'm', 'LineWidth', 1.5);
xlabel('x');
ylabel('sigma(x)');
title('Diffusion coefficient sigma(x)');
grid on;


T =min(1/abs(A)*5,100);
v = flux; sigma = diffusion;
filename ='RandFlow_v14.mat';
save(filename,'v','sigma','x','L','T');