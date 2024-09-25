% @author Joseph
% @note test sample for goldenSectionSearch.m

% Objective function for the quadratic f(x) = (x-3)^2function f = objectiveF(x)
function f = runGSS(x)
  f = (x - 3)^2; % feel free inserting your unique function here
endfunction

l = 0; u = 6; % lower and upper bonds respective
eps = 1e-5;   % tolerance level

% Call the goldenSection function with the objective function and parameters
xstar = goldenSectionSearch(@runGSS, l, u, eps);

% Display the result
printf("The estimated minimum point is at x = %.6f\n", xstar);
