function [n,prob_funct,prob_grad,prob_hd_prod,l,u,x0] = problem()
    
    % ******************************************************
    %            GENERALIZED ROSENBROCK FUNCTION
    % ******************************************************
    
    
    % PROBLEM DIMENSION
    %
    % n = problem dimension
    %
    %-------------------------------------------------------------------------------------
    n = 1000;
    %-------------------------------------------------------------------------------------
    
    
    % OBJECTIVE FUNCTION
    %
    % Required input/output arguments:
    %   x (in) is the point where the objective function will be computed
    %   f (out) is the value of objective function at x
    %   status (out) is 0 if no error occurred, >0 otherwise
    %
    %-------------------------------------------------------------------------------------
    function [f,status] = funct(x)
        
        c = 1e2;
        
        f1 = x(1:n-1).*x(1:n-1);
        f = sum(f1) + n - 1e0 - 2e0*sum(x(1:n-1)) + c*sum((x(2:n)-f1).*(x(2:n)-f1));
        
        status = 0;
    
    end
    prob_funct = @funct;
    %-------------------------------------------------------------------------------------
    
    
    % GRADIENT OF THE OBJECTIVE FUNCTION
    %
    % Required input/output arguments:
    %   x (in) is the point where the gradient of the objective function will be computed
	%   g (out) is the gradient of the objective function at x
	%   status (out) is 0 if no error occurred, >0 otherwise
    %
    %-------------------------------------------------------------------------------------
    function [g,status] = grad(x)
        
        c = 1e2;
        
        g1 = (x(1:n-1).*x(1:n-1)) - x(2:n);
        g = 2e0 * ( [-ones(n-1,1); 0e0] - c*[0e0; g1] ...
                    + [x(1:n-1); 0e0].*([ones(n-1,1); 0e0] + 2e0*c*([ones(n-1,1); 0e0].*[g1; 0e0])) );
        
        status = 0;
        
    end
    prob_grad = @grad;
    %-------------------------------------------------------------------------------------
    
    
    % HESSIAN-VECTOR PRODUCT
    %
    % Required input/output arguments:
    %   goth (in) is true if x is the same point used in the previous call, false otherwise
    %             (it can be used to store values of the Hessian matrix)
	%   x (in) is the point where the Hessian-vector product will be computed to be multiplied with d
	%   d (in) is a vector whose product with the Hessian matrix is required
    %   hd (out) is the Hessian-vector product
	%   status (out) is 0 if no error occurred, >0 otherwise
	%
	% N.B. This function can be ignored (in the sense that it can return any dummy value)
	% if Hessian-vector products are approximated (i.e., if 'hd_exact' is equal to false
	% in function 'asa_bcp')
    %
    %-------------------------------------------------------------------------------------
    function [hd,status] = hd_prod(goth,x,d)
        
        global hd1 hd2
        
        c = 1e2;
        
        if (~goth)
           hd2 = -4e0*c*x(1:n-1);
           hd1 = [2e0 + 12e0*c*x(1)*x(1); 2e0*(c+1e0)*ones(n-2,1) + 12e0*c*x(2:n-1).*x(2:n-1)] + ...
               [hd2(2:n-1); -4e0*c*x(n)];
        end
                
        hd = [hd1.*d(1:n-1) + hd2.*d(2:n); 2e0*c*d(n)] + [0e0; hd2.*d(1:n-1)];
        
        status = 0;
    
    end
    prob_hd_prod = @hd_prod;
    %-------------------------------------------------------------------------------------
    
    
    % BOUNDS
    %
    % l: lower bound
    % u: upper bound
    %
    %-------------------------------------------------------------------------------------
    l = repmat([-15e-1;5e-1],n/2,1);
    u = repmat([5e-1;2e0],n/2,1);
    if (floor(n/2) < n/2)
        l(n,1) = -15e-1;
        u(n,1) = 5e-1;
    end
    %-------------------------------------------------------------------------------------
    
    
    % STARTING POINT
    %
    % x0: starting point
    %
    %-------------------------------------------------------------------------------------
    x0 = repmat([-12e-1;1e0],floor(n/2),1);
    if (floor(n/2) < n/2)
        x0(n,1) = -12e-1;
    end
    %-------------------------------------------------------------------------------------
    
end