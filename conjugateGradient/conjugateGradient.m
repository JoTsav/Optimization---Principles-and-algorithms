% Conjugate Gradient Algorithm for Quadratic Problems
% Solves the problem: min_x (1/2 * x' * Q * x + b' * x)
% @author Joseph 	

% @param Q: Symmetric positive-definite matrix (n x n)
% @param b: Vector of size n
% @param x0: Initial guess (vector of size n)
% @param printlevel: Verbosity level (default: 0)

% @return D: Matrix of all search directions
% @return solution: Optimized solution vector

% @note run 'runCG.m' for sample results 

function [D, solution] = conjugateGradient(Q, b, x0, printlevel=0)
  n = size(x0, 1);  % Dimension of the problem
  xk = x0;          % Initial solution guess
  gk = Q * xk + b;  % Initial gradient
  dk = -gk;         % Initial search direction
  D = dk;           % Store search directions
  betak = 0;        % Initialize beta coefficient

  for k = 1:n
    if printlevel ~= 0
      fprintf('%3d\t%+10.5e\t%+10.5e\t%+10.5e\n', k, xk(1), gk(1), dk(1));
    end
    denom = dk' * Q * dk;  % Denominator of alpha
    if denom <= 0
      error('Q must be positive definite');
    end
    
    alphak = -(dk' * gk) / denom;  % Step size
    xk = xk + alphak * dk;         % Update solution
    gkp1 = Q * xk + b;             % Compute new gradient
    
    if printlevel ~= 0
      if betak == 0
        fprintf('\t%+10.5e\n', alphak);
      else
        fprintf('\t%+10.5e\t%+10.5e\n', alphak, betak);
      end
      for i = 2:n
        fprintf('\t%+10.5e\t%+10.5e\t%+10.5e\n', xk(i), gk(i), dk(i));
      end
    end
    
    betak = (gkp1' * gkp1) / (gk' * gk);  % Compute beta
    dk = -gkp1 + betak * dk;              % Update search direction
    D = [D, dk];                          % Store direction
    gk = gkp1;                            % Update gradient
  end

  if printlevel ~= 0
    fprintf('%3d\t%+10.5e\t%+10.5e\n', n+1, xk(1), gk(1));
    for i = 2:n
      fprintf('\t%+10.5e\t%+10.5e\n', xk(i), gk(i));
    end
  end
  
  solution = xk;  % returns optimal solution
end

