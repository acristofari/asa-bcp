% -------------------------------------------------------------------------
% 
% This file is part of ASA-BCP, which is a solver for bound-constrained
% optimization problems of the following form:
% 
%                                min f(x)
%                           s.t. l <= x <= u
% 
% where f(x) is a twice continuously differentiable.
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
% March 3rd, 2021
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
% Copyright 2017-2021 Andrea Cristofari, Marianna De Santis,
% Stefano Lucidi, Francesco Rinaldi.
% 
% -------------------------------------------------------------------------


function [x,f,asa_bcp_info] = asa_bcp(obj,x,l,u,opts)
    
    if (nargin < 4)
        error('at least four input arguments are required');
    end
    if (nargin > 5)
        error('at most five input arguments are required');
    end
    
    if (~isscalar(obj) || ~isstruct(obj))
        error('the first input must be a structure');
    end
    if (~isnumeric(x) || ~isreal(x) || size(x,2)>1)
        error('the second input must be a column vector');
    end
    if (~isnumeric(l) || ~isreal(l) || size(l,2)>1)
        error('the third input must be a column vector');
    end
    if (~isnumeric(u) || ~isreal(u) || size(u,2)>1)
        error('the fourth input must be a column vector');
    end
    
    n = size(x,1);
    
    if (size(l,1) ~= n)
        error('lower bound dimension must agree with problem dimension');
    end
    if (size(u,1) ~= n)
        error('upper bound dimension must agree with problem dimension');
    end
    
    % set parameters
    eps_opt = 1e-5;
    min_gd = 1e-15;
    min_norm_proj_d = 1e-9;
    min_stepsize = 1e-20;
    max_it = 1000000;
    max_n_f = 1000000;
    max_n_g = 1000000;
    max_n_hd = 1000000;
    min_f = -1e90;
    m = 100;
    z = 20;
    hd_exact = true;
    verbosity = 1;
    if (nargin == 5)
        if (~isscalar(opts) || ~isstruct(opts))
            error('the third input (which is optional) must be a structure.');
        end
        struct_field = fieldnames(opts);
        for i = 1:length(struct_field)
            switch char(struct_field(i))
                case 'eps_opt'
                    eps_opt = opts.eps_opt;
                    if (~isscalar(eps_opt) || ~isnumeric(eps_opt) || ~isreal(eps_opt) || eps_opt<=0e0)
                       error('error when calling asa_bcp: ''eps_opt'' must be greater than or equal to 0.');
                    end
                case 'min_gd'
                    min_gd = opts.min_gd;
                    if (~isscalar(min_gd) || ~isnumeric(min_gd) || ~isreal(min_gd) || min_gd<0e0)
                       error('error when calling asa_bcp: ''min_gd'' must be greater than or equal to 0.');
                    end
                case 'min_norm_proj_d'
                    min_norm_proj_d = opts.min_norm_proj_d;
                    if (~isscalar(min_norm_proj_d) || ~isnumeric(min_norm_proj_d) || ~isreal(min_norm_proj_d) || min_norm_proj_d<0e0)
                       error('error when calling asa_bcp: ''min_norm_proj_d'' must be greater than or equal to 0.');
                    end
                case 'min_stepsize'
                    min_stepsize = opts.min_stepsize;
                    if (~isscalar(min_stepsize) || ~isnumeric(min_stepsize) || ~isreal(min_stepsize) || min_stepsize<0e0)
                       error('error when calling asa_bcp: ''min_stepsize'' must be greater than or equal to 0.');
                    end
                case 'max_it'
                    max_it = opts.max_it;
                    if (~isscalar(max_it) || ~isnumeric(max_it) || ~isreal(max_it) || max_it<0e0)
                       error('error when calling asa_bcp: ''max_it'' must be greater than or equal to 0.');
                    end
                    max_it = floor(max_it);
                case 'max_n_f'
                    max_n_f = opts.max_n_f;
                    if (~isscalar(max_n_f) || ~isnumeric(max_n_f) || ~isreal(max_n_f) || max_n_f<1e0)
                       error('error when calling asa_bcp: ''max_n_f'' must be greater than or equal to 1.');
                    end
                    max_n_f = floor(max_n_f);
                case 'max_n_g'
                    max_n_g = opts.max_n_g;
                    if (~isscalar(max_n_g) || ~isnumeric(max_n_g) || ~isreal(max_n_g) || max_n_g<1e0)
                       error('error when calling asa_bcp: ''max_n_g'' must be greater than or equal to 1.');
                    end
                    max_n_g = floor(max_n_g);
                case 'max_n_hd'
                    max_n_hd = opts.max_n_hd;
                    if (~isscalar(max_n_hd) || ~isnumeric(max_n_hd) || ~isreal(max_n_hd) || max_n_hd<0e0)
                       error('error when calling asa_bcp: ''max_n_hd'' must be greater than or equal to 0.');
                    end
                    max_n_hd = floor(max_n_hd);
                case 'min_f'
                    min_f = opts.min_f;
                    if (~isscalar(min_f) || ~isnumeric(min_f) || ~isreal(min_f))
                       error('error when calling asa_bcp: ''min_f'' must be a real number.');
                    end
                case 'm'
                    m = opts.m;
                    if (~isscalar(m) || ~isnumeric(m) || ~isreal(m) || m<0e0)
                       error('error when calling asa_bcp: ''m'' must be greater than or equal to 0.');
                    end
                    m = floor(m);
                case 'z'
                    z = opts.mm;
                    if (~isscalar(z) || ~isnumeric(z) || ~isreal(z) || z<0e0)
                       error('error when calling asa_bcp: ''z'' must be greater than or equal to 0.');
                    end
                    z = floor(z);
                case 'hd_exact'
                    hd_exact = opts.hd_exact;
                    if (~isscalar(hd_exact) || ~islogical(hd_exact))
                       error('error when calling asa_bcp: ''hd_exact)'' must be a logical value.');
                    end
                case 'verbosity'
                    verbosity = opts.verbosity;
                    if (~isscalar(verbosity) || ~isnumeric(verbosity) || ~isreal(verbosity) || verbosity<0e0)
                       error('error when calling asa_bcp: ''verbosity'' must be between 0 and 2.');
                    end
                otherwise
                    error(['error when calling asa_bcp: in the fifth input argument (which is optional), ' ...
                        '''' char(struct_field(i)) ''' is not a valid field name.']);
            end
        end
    end
    
    % set objective function, gradient and Hessian-vector product
    struct_field = fieldnames(obj);
    obj_set = false;
    grad_set = false;
    hd_prod_set = false;
    for i = 1:length(struct_field)
       switch char(struct_field(i))
           case 'funct'
               if (~isa(obj.funct,'function_handle'))
                  error('error when calling asa_bcp: obj.funct must be a function handle.');
               end
               obj_set = true;
           case 'grad'
               if (~isa(obj.grad,'function_handle'))
                  error('error when calling asa_bcp: obj.grad must be a function handle.');
               end
               grad_set = true;
           case 'hd_prod'
               if (~isa(obj.hd_prod,'function_handle'))
                  error('error when calling asa_bcp: obj.hd_prod must be a function handle.');
               end
               hd_prod_set = true;
           otherwise
               error(['error when calling asa_bcp: in the first input argument, ''' char(struct_field(i)) ''' is not a valid field name.']);
       end
    end
    if (~obj_set)
        error('error when calling asa_bcp: the objective function must be specified.');
    end
    if (~grad_set)
        error('error when calling asa_bcp: the gradient of the objective must be specified.');
    end
    if (hd_exact && ~hd_prod_set)
        error(['error when calling asa_bcp: the Hessian-vector product must be specified ' ...
               '(set opts.hd_exact to false for approximating Hessian-vector products).']);
    end
    
    clear struct_field obj_set grad_set hd_prod_set
    
    l_min = -1e20; % no i-th lower bound if l(i) <= l_min
    u_max = 1e20; % no i-th upper bound if u(i) >= u_max
    
    if (verbosity > 0)
        fid = fopen('iteration_history.txt','w');
    end
    
    % identify fixed variables
    ind_true = find(l<u);
    n_true = length(ind_true);
        
    % check feasibility
    if (any(l>u))        
        x = NaN;
        f = NaN;
        sup_norm_proj_g = NaN;
        it = 0;
        n_f = 0;
        n_g = 0;
        n_hd = 0;
        it_cg_tot = 0;
        flag = -1;
        fprintf('infeasible problem\n');
        if (verbosity > 0)
            fprintf(fid,'infeasible problem\n');
        end
        return;
    end
    
    %  project the starting point onto the box
    x(l>max(l_min,x)) = l(l>max(l_min,x));
    x(u<min(u_max,x)) = u(u<min(u_max,x));
    
    % the point obtained by setting the estimated active variables to the bounds
    % is accepted in case of sufficient decrease in the objective function
    % if the distance between the two points is less than or equal to 'delta_act'
    delta_act = min(1e30,max(1e3,norm(u(ind_true(u(ind_true)<u_max & l(ind_true)>l_min)) ...
                                      -l(ind_true(u(ind_true)<u_max & l(ind_true)>l_min)))));
    
	% the unit stepsize is accepted without evaluating the objective function
    % if the distance between the two points is less than or equal to 'delta_dir'
	delta0_dir = 1e3;
    delta_dir = delta0_dir;
    beta_dir = 9e-1; % reduction factor of 'delta_dir' (must be >=0 and <1)
    
    if (verbosity > 0)
        fprintf('\nnumber of variables = %i (%i fixed)',n,n-n_true);
        fprintf(fid,'number of variables = %i (%i fixed)',n,n-n_true);
    end
    
    % first function evaluation
    [f,status] = obj.funct(x);
    if (status > 0)
        f_best = NaN;
        sup_norm_proj_g_best = NaN;
        flag = 9;
        err_obj();
        return;
    end
    n_f = 1;
    f_computed = true;
    
    x_best = x;
    f_best = f;
    f_w = f;
    w = -Inf*ones(m,1);
    
    % first gradient evaluation
    [g,status] = obj.grad(x);
    if (status > 0)
        sup_norm_proj_g_best = NaN;
        flag = 10;
        err_obj();
        return;
    end
    n_g = 1;
    
    g_best = g;
    
    gd = -Inf;
    
    n_hd = 0;
    
    % compute the sup-norm of the projected gradient
    v = x - g;
    v(l>max(l_min,v)) = l(l>max(l_min,v));
    v(u<min(u_max,v)) = u(u<min(u_max,v));
    sup_norm_proj_g = norm(x-v,Inf);
    
    sup_norm_proj_g_best = sup_norm_proj_g;
    
    % initalize 'eps_act_set'
    if (sup_norm_proj_g >= eps_opt)
        eps_act_set = min(1e-6,(norm(max(l,min(x-g,u))))^(-3));
    else
        eps_act_set = 1e-6;
    end
    
    z_nm = min(n_true,z);
    
    % initialize counters
    it = 0;
    k = 0;
    it_nm = 0;
    n_hd = 0;
    it_cg_tot = 0;
    
    act_phase = false;
    gd_exit = false;
    dir_exit = false;
    stepsize_exit = false;
    min_f_exit = (f<=min_f);
    checkpoint = true;
    is_restarted = false;
    is_first_linesearch = true;
    
    % allocate multiplier vectors
    lambda = zeros(n,1);
    mu = zeros(n,1);
    
    it_cg = 0;
    n_new_act_l = 0;
    n_new_act_u = 0;
    n_non_act = 0;
    ind_new_act_l = 0;
    ind_new_act_u = 0;
    ind_non_act = 0;
    d = zeros(n,1);
    warn_sing = false;
    warn_noposdef = false;
    warn_small = false;
    warn_conjfail = false;
    warn_maxcgit = false;
    warn_grad = false;
    is_d_neg_curv = false;
    
    flag = 0;
    
    %-----------------
    % START MAIN LOOP
    %-----------------
    
    while (~converged())
        
        if (~is_restarted)
            
            %---------------------------------------------------------------
            %     MINIMIZATION STEP OVER THE ESTIMATED ACTIVE VARIABLES
            %---------------------------------------------------------------
            
            if (~act_phase)
                % active-set estimate
                act_phase = true;
                compute_multipliers();
                estimate_active_set();
            end
            
            if (verbosity > 1)
                fprintf(['\n\n--- iteration details ---\n\n' ...
                         'number of estimated active variables not at the bounds = %-i'], ...
                        n_new_act_l+n_new_act_u);
                fprintf(fid,['\n\n--- iteration details ---\n\n' ...
                             'number of estimated active variables not at the bounds = %-i'], ...
                            n_new_act_l+n_new_act_u);
            end
            
            if (act_phase)
                
                if (~f_computed)
                    [f,status] = obj.funct(x);
                    if (status > 0)
                        flag = 9;
                        err_obj();
                        return;
                    end
                    n_f = n_f + 1;
                    f_computed = true;
                    if (verbosity > 1)
                        fprintf('\nfunction control, f = %-.5e',f);
                        fprintf(fid,'\nfunction control, f = %-.5e',f);
                    end
                end
                
                f_decreased = false;
                v = x;
                while (act_phase && ~f_decreased)
                    
                    % set the estimated active variables to the bounds
                    x([ind_new_act_l;ind_new_act_u]) = [l(ind_new_act_l);u(ind_new_act_u)];
                    
                    sq_norm_d_act = (x-v)'*(x-v);
                    
                    % compute the objective function
                    [tmp,status] = obj.funct(x);
                    if (status > 0)
                        flag = 9;
                        err_obj();
                        return;
                    end
                    n_f = n_f + 1;
                    
                    % check if the objective function is sufficiently decreased
                    if (tmp <= f-(sq_norm_d_act/(2e0*eps_act_set)))
                        f_decreased = true;
                        if (sqrt(sq_norm_d_act) <= delta_act) % new point accepted
                            delta_act = delta_act*beta_dir;
                        elseif (tmp > f_w) % point not accepted
                            if (verbosity > 1)
                                fprintf('\npoint not accepted (f = %-.5e)\n\nrestart from the best point',tmp);
                                fprintf(fid','\npoint not accepted (f = %-.5e)\n\nrestart from the best point',tmp);
                            end
                            restart();
                        end
                        if (~is_restarted) % N.B. 'f_computed' remains true
                            f = tmp;
                            % compute the gradient and the sup-norm of the projected gradient
                            [g,status] = obj.grad(x);
                            if (status > 0)
                                flag = 10;
                                err_obj();
                                return;
                            end
                            n_g = n_g + 1;
                            compute_sup_norm_proj_g();
                            if (f < f_best)
                                if (f <= min_f)
                                    min_f_exit = true;
                                else
                                    x_best = x;
                                    f_best = f;
                                    g_best = g;
                                    sup_norm_proj_g_best = sup_norm_proj_g;
                                end
                            end
                            if (verbosity > 1)
                                fprintf('\npoint accepted (f = %-.5e)',f);
                                fprintf(fid,'\npoint accepted (f = %-.5e)',f);
                            end
                        end
                    else % set x to the previous value and reduce 'eps_act_set' to carry out a new active-set estimate
                        x = v;
                        if (n_f < max_n_f)
                            eps_act_set = 1e-1*eps_act_set;
                            estimate_active_set();
                            if (verbosity > 1)
                                fprintf(['\npoint not accepted (f = %-.5e)\nreducing epsilon' ...
                                         '\nnumber of estimated active variables not at the bounds = %-i'], ...
                                        tmp,n_new_act_l+n_new_act_u);
                                fprintf(fid,['\npoint not accepted (f = %-.5e)\nreducing epsilon' ...
                                             '\nnumber of estimated active variables not at the bounds = %-i'], ...
                                            tmp,n_new_act_l+n_new_act_u);
                            end
                        else
                            act_phase = false;
                        end
                    end
                    
                end
                
            end
            
        else
            
            is_restarted = false;
            if (verbosity > 1)
                fprintf('\n\n--- iteration details ---');
                fprintf(fid,'\n\n--- iteration details ---');
            end
            
        end
        
        %---------------------------------------------------------------
        %   MINIMIZATION STEP OVER THE ESTIMATED NON-ACTIVE VARIABLES
        %---------------------------------------------------------------
        
        if (~is_restarted && sup_norm_proj_g>eps_opt && n_f<max_n_f && n_g<max_n_g && ~min_f_exit)
            
            if (act_phase)
                % active-set estimate
                act_phase = false;
                compute_multipliers();
                estimate_active_set();
            end
            
            if (verbosity > 1)
                fprintf('\n\nnumber of estimated non-active variables = %-i',n_non_act);
                fprintf(fid,'\n\nnumber of estimated non-active variables = %-i',n_non_act);
            end
            
            if (~act_phase)
                
                % compute the norm of the gradient with respect to the esimated non-active variables
                norm_g_non_act = norm(g(ind_non_act));
                
                if (norm_g_non_act > min_norm_proj_d)
                    
                    if (checkpoint)
                        % N.B. 'f_computed' = true at this point
                        if (f < f_best)
                            if (f <= min_f)
                                min_f_exit = true;
                            else
                                x_best = x;
                                f_best = f;
                                g_best = g;
                                sup_norm_proj_g_best = sup_norm_proj_g;
                            end
                        end
                        if (~min_f_exit)
                            update_w();
                            checkpoint = false;
                        end
                    end
                    
                    % function control
                    if (it_nm >= z_nm)
                        z_nm = min(z_nm+n_true,z);
                        if (~f_computed)
                            [f,status] = obj.funct(x);
                            if (status > 0)
                                flag = 9;
                                err_obj();
                                return;
                            end
                            n_f = n_f + 1;
                        end
                        if (f >= f_w)
                            if (verbosity > 1)
                                fprintf(['\nfunction control not satisfied (f = %-.5e)' ...
                                         '\n\nrestart from the best point'],f);
                                fprintf(fid,['\nfunction control not satisfied (f = %-.5e)' ...
                                             '\n\nrestart from the best point'],f);
                            end
                            restart();
                        else
                            if (verbosity > 1)
                                fprintf('\nfunction control satisfied (f = %-.5e)',f);
                                fprintf(fid,'\nfunction control satisfied (f = %-.5e)',f);
                            end
                            if (f < f_best)
                                if (f <= min_f)
                                    min_f_exit = true;
                                else
                                    x_best = x;
                                    f_best = f;
                                    g_best = g;
                                    sup_norm_proj_g_best = sup_norm_proj_g;
                                end
                            end
                            if (~min_f_exit)
                                update_w();
                                f_computed = true;
                            end
                        end
                    end
                    
                    if (~is_restarted && ~min_f_exit)
                        
                        if (hd_exact)
                            max_inner_it = min(2*n_non_act,max_n_hd-it_cg_tot);
                        else
                            max_inner_it = min(2*n_non_act,max_n_g-it_cg_tot);
                        end
                        
                        if (max_inner_it > 1)
                            
                            % compute the search direction
                            trunc_newton_dir();
                            
                            if (flag==10 || flag==11)
                                err_obj();
                                return;
                            end
                            
                            it_cg_tot = it_cg_tot + it_cg;
                            
                            if (verbosity > 1)
                                fprintf('\nnumber of inner cg iterations = %-i',it_cg);
                                fprintf(fid,'\nnumber of inner cg iterations = %-i',it_cg);
                                if (warn_sing)
                                    fprintf('\nHessian matrix probably singular');
                                    fprintf(fid,'\nHessian matrix probably singular');
                                elseif (warn_noposdef)
                                    fprintf('\nHessian matrix probably not positive definite');
                                    fprintf(fid,'\nHessian matrix probably not positive definite');
                                elseif (warn_small)
                                    fprintf('\nnorm of the inner conjugate direction too small');
                                    fprintf(fid,'\nnorm of the inner conjugate direction too small');
                                elseif (warn_conjfail)
                                    fprintf('\nconjugacy failure when computing the search direction');
                                    fprintf(fid,'\nconjugacy failure when computing the search direction');
                                elseif (warn_maxcgit)
                                    fprintf('\nconjugate gradient method not converged');
                                    fprintf(fid,'\nconjugate gradient method not converged');
                                end
                                if (warn_grad)
                                    fprintf('\nanti-gradient used as search direction');
                                    fprintf(fid,'\nanti-gradient used as search direction');
                                end
                                if (is_d_neg_curv)
                                    fprintf('\nnegative curvature direction');
                                    fprintf(fid,'\nnegative curvature direction');
                                end
                                fprintf('\ndirectional derivative = %-.5e',gd);
                                fprintf(fid,'\ndirectional derivative = %-.5e',gd);
                            end
                            
                            % check if the directional derivative is sufficiently negative
                            if (gd < -min_gd)
                                
                                % check if the line search must be performed
                                if (warn_noposdef || is_first_linesearch)
                                    ls = true;
                                else
                                    ls = false;
                                end
                                
                                % try accepting the unit stepsize or prepare for the line search
                                if (~ls)
                                    
                                    % set v = x and x = p[x + d]
                                    v = x;
                                    x = x + d;
                                    x(l>max(l_min,x)) = l(l>max(l_min,x));
                                    x(u<min(u_max,x)) = u(u<min(u_max,x));
                                    norm_proj_d = norm(x-v);
                                    
                                    if (verbosity > 1)
                                        fprintf('\nnorm of the projected direction = %-.5e',norm_proj_d);
                                        fprintf(fid,'\nnorm of the projected direction = %-.5e',norm_proj_d);
                                    end
                                    
                                    % check if the norm of the projected direction if sufficiently large
                                    if (norm_proj_d > min_norm_proj_d)
                                        if (norm_proj_d <= delta_dir) % unit stepsize accepted
                                            delta_dir = beta_dir*delta_dir;
                                            % compute the gradient and the sup-norm of the projected gradient
                                            [g,status] = obj.grad(x);
                                            if (status > 0)
                                                flag = 10;
                                                err_obj();
                                                return;
                                            else
                                                n_g = n_g + 1;
                                                compute_sup_norm_proj_g();
                                                f_computed = false;
                                                it_nm = it_nm + 1;
                                                if (verbosity > 1)
                                                    fprintf('\nstepsize = 1 accepted without computing f');
                                                    fprintf(fid,'\nstepsize = 1 accepted without computing f');
                                                end
                                            end
                                        elseif (~f_computed) % check the objective function
                                            % restore x
                                            x2 = x;
                                            x = v;
                                            v = x2;
                                            [f,status] = obj.funct(x); % evaluate f(x)
                                            if (status > 0)
                                                flag = 9;
                                                err_obj();
                                                return;
                                            end
                                            n_f = n_f + 1;
                                            if (f < f_w) % objective function decreased -> line search
                                                if (f < f_best)
                                                    if (f <= min_f)
                                                        min_f_exit = true;
                                                    else
                                                        x_best = x;
                                                        f_best = f;
                                                        g_best = g;
                                                        sup_norm_proj_g_best = sup_norm_proj_g;
                                                    end
                                                end
                                                if (~min_f_exit)
                                                    update_w();
                                                    f_computed = true;
                                                    ls = true;
                                                    % set again v = x and x = p[x + d]
                                                    v = x;
                                                    x = x2;
                                                    if (verbosity > 1)
                                                        fprintf('\nfunction control satisfied (f = %-.5e)',f);
                                                        fprintf(fid,'\nfunction control satisfied (f = %-.5e)',f);
                                                    end
                                                end
                                            else % objective function not decreased -> restart
                                                if (verbosity > 1)
                                                    fprintf(['\nfunction control not satisfied (f = %-.5e)' ...
                                                             '\n\nrestart from the best point'],f);
                                                    fprintf(fid,['\nfunction control not satisfied (f = %-.5e)' ...
                                                                 '\n\nrestart from the best point'],f);
                                                end
                                                restart();
                                            end
                                        else
                                            ls = true;
                                        end
                                    else
                                        % restore x
                                        x = v;
                                        dir_exit = true;
                                    end
                                    
                                elseif (~f_computed)
                                    
                                    [f,status] = obj.funct(x); % evaluate f(x)
                                    if (status > 0)
                                        flag = 9;
                                        err_obj();
                                        return;
                                    end
                                    n_f = n_f + 1;
                                    if (f < f_w) % objective function decreased -> line search
                                        % set v = x and x = p[x + d]
                                        f_computed = true;
                                        v = x;
                                        x = x + d;
                                        x(l>max(l_min,x)) = l(l>max(l_min,x));
                                        x(u<min(u_max,x)) = u(u<min(u_max,x));
                                        norm_proj_d = norm(x-v);
                                        if (verbosity > 1)
                                            fprintf('\nnorm of the projected direction = %-.5e',norm_proj_d);
                                            fprintf(fid,'\nnorm of the projected direction = %-.5e',norm_proj_d);
                                        end
                                        if (norm_proj_d > min_norm_proj_d)
                                            if (verbosity > 1)
                                                fprintf('\nfunction control satisfied (f = %-.5e)',f);
                                                fprintf(fid,'\nfunction control satisfied (f = %-.5e)',f);
                                            end
                                            if (f < f_best)
                                                if (f <= min_f)
                                                    x = v;
                                                    min_f_exit = true;
                                                    ls = false; % to skip the line search and exit the while loop
                                                else
                                                    x_best = x;
                                                    f_best = f;
                                                    g_best = g;
                                                    sup_norm_proj_g_best = sup_norm_proj_g;
                                                end
                                            end
                                            if (~min_f_exit)
                                                update_w();
                                            end
                                        else
                                            % restore x
                                            x = v;
                                            dir_exit = true;
                                            ls = false; % to skip the line search and exit the while loop
                                        end
                                    else % objective function not decreased -> restart
                                         % (first check if the norm of the projected direction is sufficiently large)
                                        v = x + d;
                                        v(l>max(l_min,v)) = l(l>max(l_min,v));
                                        v(u<min(u_max,v)) = u(u<min(u_max,v));
                                        norm_proj_d = norm(x-v);
                                        if (verbosity > 1)
                                            fprintf('\nnorm of the projected direction = %-.5e',norm_proj_d);
                                            fprintf(fid,'\nnorm of the projected direction = %-.5e',norm_proj_d);
                                        end
                                        if (norm_proj_d > min_norm_proj_d)
                                            if (verbosity > 1)
                                                fprintf(['\npoint not accepted (f = %-.5e)' ...
                                                         '\n\nrestart from the best point'],f);
                                                fprintf(fid,['\npoint not accepted (f = %-.5e)' ...
                                                             '\n\nrestart from the best point'],f);
                                            end
                                            restart();
                                        else
                                            dir_exit = true;
                                            f_computed = true;
                                            ls = false; % to skip the line search and exit the while loop
                                        end
                                    end
                                    
                                else
                                    
                                    % set v = x and x = p[x + d]
                                    v = x;
                                    x = x + d;
                                    x(l>max(l_min,x)) = l(l>max(l_min,x));
                                    x(u<min(u_max,x)) = u(u<min(u_max,x));
                                    norm_proj_d = norm(x-v);
                                    if (verbosity > 1)
                                        fprintf('\nnorm of the projected direction = %-.5e',norm_proj_d);
                                        fprintf(fid,'\nnorm of the projected direction = %-.5e',norm_proj_d);
                                    end
                                    if (norm_proj_d <= min_norm_proj_d)
                                        % restore x
                                        x = v;
                                        dir_exit = true;
                                        ls = false; % to skip the line search and exit the while loop
                                    end
                                    
                                end
                                
                                % line search
                                if (ls && ~is_restarted)
                                    if (verbosity >= 2)
                                        fprintf('\nline search');
                                        fprintf(fid,'\nline search');
                                    end
                                    if (is_first_linesearch)
                                        f_newton_first = f;
                                        delta_f0 = 0e0;
                                    end
                                    linesearch(); % N.B. 'f_computed' = true before and after line search
                                    if (flag == 9)
                                        err_obj();
                                        return;
                                    end
                                    if (flag ~= 5)
                                        if (~stepsize_exit) % check if the line search succeeded
                                            [g,status] = obj.grad(x);
                                            if (status > 0)
                                                flag = 10;
                                                err_obj();
                                                return;
                                            end
                                            n_g = n_g + 1;
                                            compute_sup_norm_proj_g();
                                            checkpoint = true;
                                            if (is_first_linesearch)
                                                delta_f0 = f_newton_first - f;
                                                delta_dir = delta0_dir*stepsize*norm_proj_d;
                                                is_first_linesearch = false;
                                            end
                                            if (verbosity > 1)
                                                fprintf('\nstepsize = %-.5e', stepsize);
                                                fprintf(fid,'\nstepsize = %-.5e', stepsize);
                                            end
                                        else
                                            % restore x and f(x)
                                            x = v;
                                            f = fv;
                                        end
                                    else
                                        % restore x and f(x)
                                        x = v;
                                        f = fv;
                                    end
                                end
                            else
                                gd_exit = true;
                            end
                        end
                    end
                else
                    dir_exit = true;
                end
                
            else
                
                gd = -Inf;
                norm_proj_d = Inf;
                stepsize = Inf;
                it_nm = it_nm + 1;
                
            end
            
        else
            
            gd = -Inf;
            norm_proj_d = Inf;
            stepsize = Inf;
            
        end
        
        it = it + 1;
        k = k + 1;
        
    end
    
    asa_bcp_info.sup_norm_proj_g = sup_norm_proj_g;
    asa_bcp_info.it = it;
    asa_bcp_info.inner_it = it_cg_tot;
    asa_bcp_info.n_f = n_f;
    asa_bcp_info.n_g = n_g;
    asa_bcp_info.n_hd = n_hd;
    asa_bcp_info.flag = flag;
    
    if (verbosity > 0)
        fprintf('\n\nWARNING: using ''verbosity = 0'' may be faster\n\n');
        fprintf(fid,'\n\nWARNING: using ''verbosity = 0'' may be faster\n');
        fclose(fid);
    end
    
    %-------------------------------------------------------------------------------------
    
    
    
    
    % other functions
    %-------------------------------------------------------------------------------------
    
    
    %-------------------------------------------------------------------------------------
    function conv = converged()
        if ( (sup_norm_proj_g > eps_opt) && (~gd_exit) && (~dir_exit)  && (~stepsize_exit) ...
             && (it<max_it) && (n_f<max_n_f) && (n_g<max_n_g) && (n_hd<max_n_hd) && (~min_f_exit))
            if (verbosity > 0) 
                main_prints();
            end
            conv = false;
        else 
            if (~f_computed)
                [f,status] = obj.funct(x);
                if (status > 0)
                    x = x_best;
                    f = f_best;
                    sup_norm_proj_g = sup_norm_proj_g_best;
                    flag = 9;
                    fprintf('\n\nerror when computing the objective function, algorithm stopped');
                    if (verbosity > 0)
                        fprintf(fid,'\n\nerror when computing the objective function, algorithm stopped');
                    end
                    conv = true;
                    return;
                end
                n_f = n_f + 1;
                f_computed = true;
                if (f <= min_f)
                    flag = 8;
                    if (verbosity > 0)
                        fprintf('\n\nobjective value below the minimum threshold, algorithm stopped');
                        fprintf(fid,'\n\nobjective value below the minimum threshold, algorithm stopped');
                    end
                    conv = true;
                    return;
                end
            end
            if (sup_norm_proj_g <= eps_opt)
                if (verbosity > 0)
                    main_prints();
                    fprintf(['\n\n=================================================================================' ...
                             '\noptimality condition satisfied: sup-norm of the projected gradient <= %-.5e' ...
                             '\n================================================================================='],eps_opt);
                    fprintf(fid,['\n\n=================================================================================' ...
                                 '\noptimality condition satisfied: sup-norm of the projected gradient <= %-.5e' ...
                                 '\n================================================================================='],eps_opt);
                end
                if (f <= f_best)
                    flag = 0;
                    conv = true;
                else 
                    conv = check_stop(true);
                end
            elseif (gd_exit)
                if (verbosity > 0)
                    fprintf('\n\ndirectional derivative not sufficiently negative');
                    fprintf(fid,'\n\ndirectional derivative not sufficiently negative');
                end
                % active-set estimate
                act_phase = true;
                compute_multipliers();
                estimate_active_set();                
                if (act_phase)
                    gd_exit = false;
                    conv = check_stop(false);
                else
                    if (f <= f_best)
                        it = it - 1;
                        flag = 1;
                        if (verbosity > 0)
                            fprintf(', algorithm stopped');
                            fprintf(fid,', algorithm stopped');
                        end
                        conv = true;
                    else
                        conv = check_stop(true);
                    end
                end
            elseif (dir_exit)
                if (verbosity > 0)
                    fprintf('\n\nnorm of the projected direction too small');
                    fprintf(fid,'\n\nnorm of the projected direction too small');
                end
                % active-set estimate
                act_phase = true;
                compute_multipliers();
                estimate_active_set();                
                if (act_phase)
                    dir_exit = false;
                    conv = check_stop(false);
                else
                    if (f <= f_best)
                        it = it - 1;
                        flag = 2;
                        if (verbosity > 0)
                            fprintf(', algorithm stopped');
                            fprintf(fid,', algorithm stopped');
                        end
                        conv = true;
                    else
                        conv = check_stop(true);
                    end
                end
            elseif (stepsize_exit)
                if (verbosity > 0)
                    fprintf('\n\nstepsize too small');
                    fprintf(fid,'\n\nstepsize too small');
                end
                % active-set estimate
                act_phase = true;
                compute_multipliers();
                estimate_active_set();                
                if (act_phase)
                    stepsize_exit = false;
                    conv = check_stop(false);
                else
                    if (f <= f_best)
                        it = it - 1;
                        flag = 3;
                        if (verbosity > 0)
                            fprintf(', algorithm stopped');
                            fprintf(fid,', algorithm stopped');
                        end
                        conv = true;
                    else
                        conv = check_stop(true);
                    end
                end
            elseif (min_f_exit)
                flag = 8;
                if (verbosity > 0)
                    fprintf('\n\nobjective value below the minimum threshold, algorithm stopped');
                    fprintf(fid,'\n\nobjective value below the minimum threshold, algorithm stopped');
                end
                conv = true;
            else
                if (f > f_best)
                    x = x_best;
                    f = f_best;
                    sup_norm_proj_g = sup_norm_proj_g_best;
                end
                if (it >= max_it)
                    flag = 4;
                    if (verbosity > 0)
                        main_prints();
                        fprintf('\n\ntoo many iterations, algorithm stopped');
                        fprintf(fid,'\n\ntoo many iterations, algorithm stopped');
                    end
                elseif (n_f >= max_n_f)
                    flag = 5;
                    if (verbosity > 0)
                        main_prints();
                        fprintf('\n\ntoo many function evaluations, algorithm stopped');
                        fprintf(fid,'\n\ntoo many function evaluations, algorithm stopped');
                    end
                elseif (n_g >= max_n_g)
                    flag = 6;
                    if (verbosity > 0)
                        main_prints();
                        fprintf('\n\ntoo many gradient evaluations, algorithm stopped');
                        fprintf(fid,'\n\ntoo many gradient evaluations, algorithm stopped');
                    end
                else % n_hd >= max_n_hd
                    flag = 7;
                    if (verbosity > 0)
                        main_prints();
                        fprintf('\n\ntoo many Hessian-vector products, algorithm stopped');
                        fprintf(fid,'\n\ntoo many Hessian-vector products, algorithm stopped');
                    end
                end
                conv = true;
            end
        end  
    end
    %-------------------------------------------------------------------------------------
    
    
    %-------------------------------------------------------------------------------------
    function stp = check_stop(rst) % it returns true if the algorithm must stop, false otherwise
        if ( (it<max_it) && (n_f<max_n_f) && (n_g<max_n_g) && (n_hd<max_n_hd) )
            if (rst)
                restart();
                if (verbosity > 0)
                    if (verbosity > 1)
                        fprintf('\n\nrestart from the best point');
                        fprintf(fid,'\n\nrestart from the best point');
                    end
                    main_prints();
                end
            end
            stp = false;
        else
            if (rst || f>f_best)
                x = x_best;
                f = f_best;
                sup_norm_proj_g = sup_norm_proj_g_best;
            end
            if (it >= max_it)
                flag = 4;
                if (verbosity > 0)
                    fprintf('\ntoo many iterations, algorithm stopped');
                    fprintf(fid,'\ntoo many iterations, algorithm stopped');
                end
            elseif (n_f >= max_n_f)
                flag = 5;
                if (verbosity > 0)
                    fprintf('\ntoo many function evaluations, algorithm stopped');
                    fprintf(fid,'\ntoo many function evaluations, algorithm stopped');
                end
            elseif (n_g >= max_n_g)
                flag = 6;
                if (verbosity > 0)
                    fprintf('\ntoo many gradient evaluations, algorithm stopped');
                    fprintf(fid,'\ntoo many gradient evaluations, algorithm stopped');
                end
            else % n_hd >= max_n_hd
                flag = 7;
                if (verbosity > 0)
                    fprintf('\ntoo many Hessian-vector products, algorithm stopped');
                    fprintf(fid,'\ntoo many Hessian-vector products, algorithm stopped');
                end
            end
            stp = true;
        end
    end
    %-------------------------------------------------------------------------------------
    
    
    %-------------------------------------------------------------------------------------
    function main_prints()
        fprintf(['\n\n--------------------------------------------------------\n\n' ...
                 'iteration %-i\n\n' ...
                 'best f = %-.5e\n' ...
                 'sup-norm of the projected gradient at the current point = %-.5e\n' ...
                 'number of function evaluations = %-i\n' ...
                 'number of gradient evaluations = %-i\n' ...
                 'number of Hessian-vector products = %-i\n' ...
                 'number of inner cg iterations = %-i'], ...
                it,min(f_best,f),sup_norm_proj_g,n_f,n_g,n_hd,it_cg_tot);
        fprintf(fid,['\n\n--------------------------------------------------------\n\n' ...
                     'iteration %-i\n\n' ...
                     'best f = %-.5e\n' ...
                     'sup-norm of the projected gradient at the current point = %-.5e\n' ...
                     'number of function evaluations = %-i\n' ...
                     'number of gradient evaluations = %-i\n' ...
                     'number of Hessian-vector products = %-i\n' ...
                     'number of inner cg iterations = %-i'], ...
                    it,min(f_best,f),sup_norm_proj_g,n_f,n_g,n_hd,it_cg_tot);
    end
    %-------------------------------------------------------------------------------------
    
    
    %-------------------------------------------------------------------------------------
    function compute_sup_norm_proj_g()
        proj_g = x - g;
        proj_g(l>max(l_min,proj_g)) = l(l>max(l_min,proj_g));
        proj_g(u<min(u_max,proj_g)) = u(u<min(u_max,proj_g));
        sup_norm_proj_g = norm(x-proj_g,Inf);
    end
    %-------------------------------------------------------------------------------------
     
    
    %-------------------------------------------------------------------------------------
    function compute_multipliers()
        mult1 = (u-x).*(u-x);
        mult2 = (x-l).*(x-l);
        lambda = g.*mult1./(mult1+mult2);
        mu = -g.*mult2./(mult1+mult2);
    end
    %-------------------------------------------------------------------------------------
    
    
    %-------------------------------------------------------------------------------------
    function estimate_active_set()
        est = true;
        count = 0; % to avoid cycling
        while (est && count<2)
            if (act_phase)
                ind_new_act_l = ind_true((x(ind_true)<=l(ind_true)+eps_act_set*lambda(ind_true)) & ...
                                         (g(ind_true)>0e0) & (x(ind_true)>l(ind_true)) & (l(ind_true)>l_min)) ;
                ind_new_act_u = ind_true((x(ind_true)>=u(ind_true)-eps_act_set*mu(ind_true)) & ...
                                         (g(ind_true)<0e0) & (x(ind_true)<u(ind_true)) & (u(ind_true)<u_max));
                n_new_act_l = length(ind_new_act_l);
                n_new_act_u = length(ind_new_act_u);
                act_phase = (n_new_act_l+n_new_act_u>0);
                est = ~act_phase;
                count = count + 1;
            else
                ind_non_act = ind_true(~((x(ind_true)<=l(ind_true)+eps_act_set*lambda(ind_true)) & ...
                                         (g(ind_true)>0e0) & (l(ind_true)>l_min)) & ...
                                       ~((x(ind_true)>=u(ind_true)-eps_act_set*mu(ind_true)) & ...
                                         (g(ind_true)<0e0) & (u(ind_true)<u_max)));
                n_non_act = length(ind_non_act);
                act_phase = (n_non_act==0);
                est = act_phase;
                count = count + 1;
            end
        end
    end
    %-------------------------------------------------------------------------------------
    
    
    %-------------------------------------------------------------------------------------
    function update_w() % update the vector of reference values
        w = [f; w(1:m-1)];
        f_w = max(w);
        it_nm = 0;
    end
    %-------------------------------------------------------------------------------------
    
    
    %-------------------------------------------------------------------------------------
    function restart()
        x = x_best;
        f = f_best;
        g = g_best;
        sup_norm_proj_g = sup_norm_proj_g_best;
        f_w = f_best;
        m = round((m+4)/5);
        z_nm = min(n_true,z);
        w = f_w*ones(m,1);
        k = 0;
        it_nm = 0;
        delta0_dir = delta0_dir*1e-1;
        f_computed = true;
        checkpoint = false;
        is_restarted = true;
        is_first_linesearch = true;
        act_phase = true;
    end
    %-------------------------------------------------------------------------------------
    
    
    %-------------------------------------------------------------------------------------
    function trunc_newton_dir()
        eps_cg_dir = min(1e0,min_norm_proj_d^(3e0/2e0)); % 'eps_cg_dir' should be <= 'min_norm_proj_d',
                                                         % so that the norm of the first conjugate
                                                         % direction is sufficiently large
        eps_cg_curv = 1e-9;
        
        % In this function, 'p' are the (normalized) conjugate directions in the subspace
        % of the estimated non-active variables and 'd' is the resulting search direction
        
        warn_grad = false;
        warn_small = false;
        warn_noposdef = false;
        warn_conjfail = false;
        warn_maxcgit = false;
        is_d_neg_curv = false;
        it_cg = 1;
        
        g_non_act = g(ind_non_act);
        
        if (hd_exact)
            p = zeros(n,1);
        else
            v = x;
        end
        
        % compute the first 'p' and normalize it
        p_non_act = -g_non_act/norm_g_non_act;        
        
        % compute H(x)*p
        if (hd_exact)
            p(ind_non_act) = p_non_act;
            [hp,status_tn] = obj.hd_prod(false,x,p);
            if (status_tn > 0)
                flag = 11;
                return;
            end
            n_hd = n_hd + 1;
        else
            approximate_hd();
            if (status_tn > 0)
                flag = 10;
                return;
            end
        end
        hp_non_act = hp(ind_non_act);
        
        % compute curvature
        curv = p_non_act'*hp_non_act;
        
        % check curvature
        if (curv > eps_cg_curv)
            
            % update 'd'
            d_non_act = -g_non_act/curv;
            
            % first truncated-Newton termination test
            if (curv <= 1e0)
                eps_cg_tr = min(1e0,1e-1+exp(-1e-3*double(k)));
            else
                eps_cg_tr = min(1e0,1e-1+exp(-1e-3*double(k)))*min(1e0,norm_g_non_act/curv);
            end
            if (eps_cg_tr < 1e0)
                
                gd = -norm_g_non_act*norm_g_non_act/curv;
                q = 5e-1*gd;
                
                hd_non_act = hp_non_act*(norm_g_non_act/curv);
                g_q = g_non_act + hd_non_act;
                
                % check if the maximum number of inner iterations has been reached
                if (max_inner_it > 1)
                    
                    % update 'p'
                    beta = (g_q'*hp_non_act)/curv;
                    p_non_act = -g_q + beta*p_non_act;
                    
                    gp = g_non_act'*p_non_act;
                    gqp = g_q'*p_non_act;
                    
                    % check conjugacy
                    if ((abs(gqp-gp)>abs(gqp)+1e-6) || ((gp*gqp<0e0) && (abs(gqp-gp)>1e-6*(abs(gqp)+1e-9))))
                        warn_conjfail = true;
                    end
                    
                    trunc_exit = false;
                    
                    % repeat for iteration 2,3,...
                    while (~trunc_exit && ~warn_conjfail && ~warn_noposdef && ~warn_small && ~warn_maxcgit)
                        
                        it_cg = it_cg + 1;
                        
                        norm_p = norm(p_non_act);
                        
                        % check if the norm of 'p' is sufficiently large
                        if (norm_p > eps_cg_dir)
                            
                            % normalize 'p'
                            p_non_act = p_non_act/norm_p;
                            
                            % compute H(x)*p
                            if (hd_exact)
                                p(ind_non_act) = p_non_act;
                                [hp,status_tn] = obj.hd_prod(true,x,p);
                                if (status_tn > 0)
                                    flag = 11;
                                    return;
                                end
                                n_hd = n_hd + 1;
                            else
                                status_tn = approximate_hd();
                                if (status_tn > 0)
                                    flag = 10;
                                    return;
                                end
                            end
                            hp_non_act = hp(ind_non_act);
                            
                            % compute curvature
                            curv = p_non_act'*hp_non_act;
                            
                            % check curvature
                            if (curv > eps_cg_curv)
                                
                                % update 'd'
                                alpha = -(p_non_act'*g_q)/curv;
                                d_non_act = d_non_act + alpha*p_non_act;
                                
                                g_q = g_q + alpha*hp_non_act;
                                hd_non_act = hd_non_act + alpha*hp_non_act;
                                
                                % truncated-Newton termination test
                                norm_d = norm(d_non_act);
                                if (norm_d >= norm_g_non_act)
                                    eps_cg_tr = min(1e0,1e-1+exp(-1e-3*double(k)));
                                else
                                    eps_cg_tr = min(1e0,1e-1+exp(-1e-3*double(k)))*min(1e0,norm_d);
                                end
                                q_old = q;
                                gd_old = gd;
                                gd = g_non_act'*d_non_act;
                                q = ((g_q'*d_non_act)+gd)/2e0;
                                if (~( abs((q-q_old-3e0*(gd-gd_old)/2e0) / (q-3e0*gd/2e0))*it_cg <= eps_cg_tr ))
                                    
                                    % check if the maximum number of inner iterations has been reached
                                    if (it_cg < max_inner_it)
                                        
                                        % update 'p'
                                        beta = (g_q'*hp_non_act)/curv;
                                        p_non_act = -g_q + beta*p_non_act;
                                        
                                        gp = g_non_act'*p_non_act;
                                        gqp = g_q'*p_non_act;
                                        
                                        % check conjugacy
                                        if ((abs(gqp-gp)>abs(gqp)+1e-6) || ((gp*gqp<0e0) && (abs(gqp-gp)>1e-6*(abs(gqp)+1e-9))))
                                            warn_conjfail = true;
                                        end
                                        
                                    else
                                        warn_maxcgit = true;
                                    end
                                else
                                    trunc_exit = true;
                                end
                                
                            else
                                
                                % Hessian matrix not (sufficiently) positive definite
                                warn_noposdef = true;
                                if (curv < -eps_cg_curv)
                                	alpha = -(p_non_act'*g_q)/curv;
                                    d_non_act = d_non_act - alpha*p_non_act;
                                    
                                    hd_non_act = hd_non_act - alpha*hp_non_act;
                                    
                                    gd = g_non_act'*d_non_act;
                                    curv = d_non_act'*hd_non_act;
                                    sq_norm_d = d_non_act'*d_non_act;
                                    
                                    if (curv < eps_cg_curv*sq_norm_d)
                                        is_d_neg_curv = (curv < -eps_cg_curv*sq_norm_d);
                                        warn_sing = ~is_d_neg_curv;
                                    end
                                else
                                    warn_sing = true;
                                end
                                
                            end
                        else
                            warn_small = true;
                        end
                    end
                else
                    warn_maxcgit = true;
                end
            else
                gd = -norm_g_non_act*norm_g_non_act/curv;
            end
            
        else
            
            % Hessian matrix not (sufficiently) positive definite
            warn_noposdef = true;
            
            if (curv < -eps_cg_curv)
                d_non_act = g_non_act/curv;
                gd = norm_g_non_act*norm_g_non_act/curv;
                is_d_neg_curv = true;
            else
                d_non_act = -g_non_act;
                gd = norm_g_non_act*norm_g_non_act;
                warn_sing = true;
                warn_grad = true;
            end
                 
        end
        
        d = zeros(n,1);
        d(ind_non_act) = d_non_act;
        
        function approximate_hd()
            eps_approx = 1e-6;
            v(ind_non_act) = x(ind_non_act) + eps_approx*p_non_act;
            [g1,status_tn] = obj.grad(v);
            if (status_tn > 0)
                return;
            end
            n_g = n_g + 1;
            hp = (g1-g)/eps_approx;
        end
        
    end
    %-------------------------------------------------------------------------------------
    
    
    %-------------------------------------------------------------------------------------
    function linesearch()
        gamma = 1e-6;
        delta = 5e-1;
        f_dec = 1e6;
        
        % N.B. 'v' is the current point, 'x' is the point obtained by using
        % the search direction with unit stepsize (and projecting onto the box)
        % and 'f' is f(v)
        
        red_param_computed = false;
        
        f_v = f;
        [f,status] = obj.funct(x);
        if (status > 0)
            flag = 9;
            return;
        end
        n_f = n_f + 1;
        
        if ((f-f_newton_first>=f_dec*delta_f0) || warn_noposdef)
            f_ref = f_v;
        else
            f_ref = f_w;
        end
        
        stepsize = 1e0;
        
        while (f > f_ref+gamma*stepsize*gd)
            
            if (n_f >= max_n_f)
                flag = 5;
                break;
            end
            
            % update the stepsize
            if ((f-f_v)/min(1e12,max(1e-12,-gamma*stepsize*gd)) < 1e10)
                stepsize = delta*stepsize;
            else
                if (~red_param_computed)
                    a1 = max(1e0,norm(v))/max(1e-15,norm_proj_d);
                    red_param_computed = true;
                end
                a1 = min(a1,stepsize*delta*delta);
                aa = 2e-12*stepsize;
                if (max(aa,a1) > min_stepsize)
                    stepsize = max(aa,a1);
                else
                    stepsize = delta*stepsize;
                end
            end
            
            if (stepsize <= min_stepsize)
                stepsize_exit = true;
                break;
            end
            
            % compute a new point
            x = v + stepsize*d;
            x(l>max(l_min,x)) = l(l>max(l_min,x));
            x(u<min(u_max,x)) = u(u<min(u_max,x));
            [f,status] = obj.funct(x);
            if (status > 0)
                flag = 9;
                return;
            end
            n_f = n_f + 1;
            
        end
        
    end
    %-------------------------------------------------------------------------------------
    
    
    %-------------------------------------------------------------------------------------
    function err_obj()
        x = x_best;
        f = f_best;
        sup_norm_proj_g = sup_norm_proj_g_best;
        asa_bcp_info.sup_norm_proj_g = sup_norm_proj_g;
        asa_bcp_info.it = it;
        asa_bcp_info.inner_it = it_cg_tot;
        asa_bcp_info.n_f = n_f;
        asa_bcp_info.n_g = n_g;
        asa_bcp_info.n_hd = n_hd;
        asa_bcp_info.flag = flag;
        if (flag == 9)
            fprintf('\n\nerror when computing the objective function, algorithm stopped\n\n');
            if (verbosity > 0)
                fprintf('WARNING: using ''verbosity = 0'' may be faster\n\n');
                fprintf(fid,['\n\nerror when computing the objective function, algorithm stopped\n\n' ...
                             'WARNING: using ''verbosity = 0'' may be faster\n']);
                fclose(fid);
            end
        elseif (flag == 10)
            fprintf('\n\nerror when computing the gradient of the objective function, algorithm stopped\n\n');
            if (verbosity > 0)
                fprintf('WARNING: using ''verbosity = 0'' may be faster\n\n');
                fprintf(fid,['\n\nerror when computing the gradient of the objective function, algorithm stopped\n\n' ...
                             'WARNING: using ''verbosity = 0'' may be faster\n']);
                fclose(fid);
            end
        else % flag == 11
            fprintf('\n\nerror when computing the Hessian-vector product, algorithm stopped\n\n');
            if (verbosity > 0)
                fprintf('WARNING: using ''verbosity = 0'' may be faster\n\n');
                fprintf(fid,['\n\nerror when computing the Hessian-vector product, algorithm stopped\n\n' ...
                             'WARNING: using ''verbosity = 0'' may be faster\n']);
                fclose(fid);
            end
        end
    end
    %-------------------------------------------------------------------------------------

end