

% Newton Line search
% JoTshava

function [solution, iteres, niter] = newtonLineSearch(obj, x0, eps, printlevel = 1, maxiter = 100)
  % Newton's method with line search and modified Cholesky for optimization
  % Inputs:
  % obj: Objective function handle (returns function value, gradient, Hessian)
  % x0: Initial guess for solution
  % eps: Tolerance for stopping criterion
  % printlevel: Controls output verbosity (optional)
  % maxiter: Maximum number of iterations (optional)
  % Initialize variables
  iteres = zeros(maxiter + 1, 4);  % Store iteration results
  xk = x0;  % Initial point
  [f, g, H] = feval(obj, xk);  % Function, gradient, and Hessian at x0
  iteres(1, :) = [xk' f norm(g)];  % Store initial values
  k = 0;  % Iteration counter
  
  % Print initial details if required
  if printlevel
    printf("f\t\t||g||\t\talpha\t\ttau\n");
    printf("%e\t%e\n", f, norm(g));
  endif
  % Line search and Cholesky parameters
  alpha0 = 1.0;
  beta1 = 1.0e-4;
  beta2 = 0.99;
  lambda = 2.0;
  % Main loop: continues until the gradient norm is small or max iterations reached
  while norm(g) > eps && k < maxiter
    % Modified Cholesky decomposition for direction d
    [L, tau] = modifiedCholesky(H);
    d = -L' \ (L \ g);  % Descent direction

    % Perform line search for optimal step size
    alpha = lineSearch(obj, xk, d, alpha0, beta1, beta2, lambda, printlevel);
    xk = xk + alpha * d;  % Update xk based on step size and direction
    
    % Evaluate new function, gradient, and Hessian
    [f, g, H] = feval(obj, xk);
    % Print iteration details if required
    if printlevel
      printf("%e\t%e\t%f\t%e\n", f, norm(g), alpha, tau);
    endif
    % Store iteration results
    k = k + 1;
    iteres(k + 1, :) = [xk' f norm(g)];
  endwhile
  % Return solution and iteration results
  solution = xk;
  niter = k;
endfunction
