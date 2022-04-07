% -------------------------------------------------------------------------
% 
% This file is part of ASA-BCP, which is a solver for bound-constrained
% optimization problems of the following form:
% 
%                                min f(x)
%                           s.t. l <= x <= u
% 
% with given vectors l, u and where f(x) is a twice continuously
% differentiable function.
% 
% -------------------------------------------------------------------------
% 
% Reference paper:
% 
% A. Cristofari, M. De Santis, S. Lucidi, F. Rinaldi (2017). A Two-Stage
% Active-Set Algorithm for Bound-Constrained Optimization. Journal of
% Optimization Theory and Applications, 172(2), 369-401.
% 
% -------------------------------------------------------------------------
% 
% Authors:
% Andrea Cristofari (e-mail: andrea.cristofari@unipd.it)
% Marianna De Santis (e-mail: mdesantis@diag.uniroma1.it)
% Stefano Lucidi (e-mail: lucidi@diag.uniroma1.it)
% Francesco Rinaldi (e-mail: rinaldi@math.unipd.it)
% 
% Last update of this file:
% April 7th, 2022
% 
% Licensing:
% This file is part of ASA-BCP.
% ASA-BCP is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% ASA-BCP is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
% GNU General Public License for more details.
% You should have received a copy of the GNU General Public License
% along with ASA-BCP. If not, see <http://www.gnu.org/licenses/>.
% 
% Copyright 2017-2022 Andrea Cristofari, Marianna De Santis,
% Stefano Lucidi, Francesco Rinaldi.
% 
% -------------------------------------------------------------------------


clear all, clc;

% In this file, it is shown how to call ASA-BCP to solve a user-defined problem.

% (1) Load the problem
[n,obj.funct,obj.grad,obj.hd_prod,l,u,x0] = problem();

% (2) Call ASA-BCP
[x,f,asa_bcp_info] = asa_bcp(obj,x0,l,u);

%--------------------------------------------------------------------------
% *** EXAMPLE OF HOW TO CHANGE ASA-BCP PARAMETERS ***
% (see the file 'usage.txt' to know which parameters can be changed and
% their default values)
%
% Instead of calling ASA-BCP by the above instruction, do the following:
%
% - create a structure having as field names the names of the parameters
%   to be changed and assign them new values, e.g.,
%
%     opts.verbosity = 0;
%
% - pass the structure to ASA-BCP as fifth input argument, e.g.,
%
%     [x,f,asa_bcp_info] = asa_bcp(obj,x0,l,u,opts);
%--------------------------------------------------------------------------

% write some statistics to the screen
if (flag >= 0)
    fprintf(['************************************************' ...
         '\n\nAlgorithm: ASA-BCP' ...
         '\n\nnumber of variables = %-i' ...
         '\n\nf = %-.5e' ...
         '\n\nsup-norm of the projected gradient = %-.5e' ...
         '\nnumber of iterations = %-i' ...
         '\nnumber of function evaluations = %-i' ...
         '\nnumber of gradient evaluations = %-i' ...
         '\nnumber of Hessian-vector products = %-i' ...
         '\nnumber of inner cg iterations = %-i' ...
         '\nexit flag = %-i' ...
         '\n\n************************************************\n'], ...
        n,f,asa_bcp_info.sup_norm_proj_g,asa_bcp_info.it,asa_bcp_info.n_f, ...
        asa_bcp_info.n_g,asa_bcp_info.n_hd,asa_bcp_info.inner_it,asa_bcp_info.flag);
else
    fprintf('%s\n','infeasible problem');
end