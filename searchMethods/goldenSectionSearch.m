
% Goldden Section Search
% @author Joseph 	

% Function to perform Golden Section Search to find the minimum of a given function
% Inputs:
% obj: Objective function handle
% l: Lower bound of the search interval
% u: Upper bound of the search interval
% eps: Tolerance level for termination

% @note run 'runGS.m' for sample results 

function xstar = goldenSectionSearch(obj, l, u, eps)
  

  % Golden ratio constant
  rho = (3.0 - sqrt(5.0)) / 2.0;
  
  % Initial internal points
  alpha1 = l + rho * (u - l);
  alpha2 = u - rho * (u - l);
  
  % Evaluate the objective function at the internal points
  h1 =feval(obj,alpha1);
  h2 =feval(obj,alpha2);
  
  k = 1; % Iteration counter

  % Main loop to iteratively reduce the search interval
  while ((u - l) > eps)
    % Display iteration details
    printf("%d\t%.6g\t%.6g\t%.6g\t%.6g\t%.6g\t%.6g\n", k, l, alpha1, alpha2, u, h1, h2);
    
    % Check if both function evaluations are equal
    if (h1 == h2)
      % Reduce both the upper and lower bounds
      l = alpha1;
      u = alpha2;
    elseif (h1 > h2)
      % Shift the lower bound and update internal points
      l = alpha1;
      alpha1 = alpha2;
      h1 = h2;
    else
      % Shift the upper bound and update internal points
      u = alpha2;
      alpha2 = alpha1;
      h2 = h1;
    endif
    % Recalculate internal points and objective function values
    alpha1 = l + rho * (u - l);
    alpha2 = u - rho * (u - l);
    h1 = feval(obj, alpha1);
    h2 = feval(obj, alpha2);
    
    % Increment iteration counter
    k = k + 1;
  endwhile
  %Final iteration output
  printf("%d\t%.6g\t%.6g\t%.6g\t%.6g\t%.6g\t%.6g\n", k, l, alpha1, alpha2, u, h1, h2);
  
  %Return the midpoint as the best estimate for the minimum
  xstar = (l + u) / 2.0;
endfunction

