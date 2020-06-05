% -------------------------------------------------------------------------
%
% This file is part of ASA-BCP, which is a solver for bound-constrained
% optimization problems of the following form:
%
%                                 min f(x)
%                           s.t. l <= x <= u
%
% with f(x) twice continuously differentiable.
%
% This is a driver for running ASA-BCP on user-defined problems.
% See the file 'README.txt' to know how to compile and run the program.
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
% June 5th, 2020
%
% Copyright 2017-2020 Andrea Cristofari, Marianna De Santis,
% Stefano Lucidi, Francesco Rinaldi.
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
% -------------------------------------------------------------------------

clear all, clc

% load the problem
[n,obj.funct,obj.grad,obj.hd_prod,l,u,x0] = problem();

% call the solver
t_start = tic;
[x,f,asa_bcp_info] = asa_bcp(obj,x0,l,u);
t_tot = toc(t_start);

%-------------------------------------------------------------------------------------------
% *** EXAMPLE OF HOW TO CHANGE ASA-BCP PARAMETERS ***
% (see the description of asa_bcp in file 'asa_bcp.m' to know which 
% parameters can be changed and their default values)
%
% Instead of calling the solver by the above instruction
% '[x,f,asa_bcp_info] = asa_bcp(obj,x0,l,u);', do the following:
%
% (1) create a structure having as field names the names of the parameters
%     to be changed and assign them new values, for instance:
%
%       opts.verbosity = 0;
%
% (2)  pass the structure to 'asa_bcp' as fifth input argument, for instance:
%
%       [x,f,asa_bcp_info] = asa_bcp(obj,x0,l,u,opts);
%-------------------------------------------------------------------------------------------

% write statistics to the screen and to file 'statistics.txt'
fid = fopen('statistics.txt','w');
if (flag < 0)
    fprintf(fid,'%s\n','infeasible problem');
    fclose(fid);
	return;
end
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
         '\nelapsed time (s) = %-.4e' ...
         '\n\n************************************************\n'], ...
        n,f,asa_bcp_info.sup_norm_proj_g,asa_bcp_info.it,asa_bcp_info.n_f, ...
        asa_bcp_info.n_g,asa_bcp_info.n_hd,asa_bcp_info.inner_it,asa_bcp_info.flag,max(t_tot,0e0));
fprintf(fid,['************************************************' ...
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
            '\nelapsed time (s) = %-.4e' ...
            '\n\n************************************************\n'], ...
           n,f,asa_bcp_info.sup_norm_proj_g,asa_bcp_info.it,asa_bcp_info.n_f, ...
           asa_bcp_info.n_g,asa_bcp_info.n_hd,asa_bcp_info.inner_it,asa_bcp_info.flag,max(t_tot,0e0));
fclose(fid);

% write the solution found by the algorithm to file 'opt_sol.txt'
fid = fopen('opt_sol.txt','w');
fprintf(fid,'%-.5e\n',x);
fclose(fid);