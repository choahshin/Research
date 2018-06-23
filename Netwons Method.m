%% Newton's method
% Nonlinear system solver f(x) = 0

% define a system
f = @(x) [x(1).^2 + x(2).^2 - 1; x(1) + x(2) - 1];

% Jacobian matrix
Df = @(x) [2.*x(1), 2.*x(2); 1, 1];

% choose initial guess
x0 = [.5;-1];
%x0 = [.5;1];

% tolerance
Err = 1e-5; % absolute error tolerance
Era = 1e-10; % relative error tolerance
n = 0;      % counting parameter
M = 100;    % maximum number of iteration allowed

x = x0;

% iterations
while norm(f(x),2) > Err.*norm(f(x0),2) + Era || n>M
    n = n+1; % counting parameter
    dx = -Df(x)\f(x);
    x = x + dx; % new guess
    f(x);
end
