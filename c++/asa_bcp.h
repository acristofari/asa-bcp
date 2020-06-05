#ifndef __ASA_BCP_H_INCLUDED__
#define __ASA_BCP_H_INCLUDED__
#include <vector>
#include <fstream>

class Problem;

struct asa_bcp_options {

    // *** DO NOT CHANGE DEFINITIONS IN THIS STRUCTURE ***
    //
    // to change ASA-BCP parameters, in the main function create an
    // object of structure type asa_bcp_options, assign new values to
    // (some of) its members and pass the address of the structure object
    // as third input argument when calling the Asa_bcp constructor
    // (see file 'main.cpp' for an example).


    // ====================================
    // DEFAULT VALUES OF ASA-BCP PARAMETERS
    // ====================================

    // PARAMETERS FOR TERMINATION (see the description of ASA_BCP in file 'asa_bcp.cpp')
    double eps_opt = 1e-5;
    double min_gd = 1e-15;
    double min_norm_proj_d = 1e-9;
    double min_stepsize = 1e-20;
    double min_f = -1e90;
    int max_it = 1000000;
    int max_n_f = 1000000;
    int max_n_g = 1000000;
    int max_n_hd = 1000000;

    // OTHER ALGORITHM PARAMETERS (see the description of ASA_BCP in file 'asa_bcp.cpp')
    int m = 100;
    int z = 20;
    bool hd_exact = true;
    short int verbosity = 1;

};

class Asa_bcp {

private:
    Problem *prob;

    unsigned int max_it,max_n_f,max_n_g,max_n_hd,m,z;
    unsigned int n,n_true,z_nm,it,k,it_cg,it_cg_tot,it_nm,n_non_act,n_new_act_l,n_new_act_u,j,h,max_inner_it,n_f,n_g,n_hd;
    short int flag;
    unsigned short int verbosity,status;
    std::vector<unsigned int> ind_new_act,ind_non_act;
    double eps_opt,min_gd,min_norm_proj_d,min_stepsize,min_f,delta0_dir;
    double f,f_best,sup_norm_proj_g_best,f_w,fv,stepsize,eps_act_set,norm_g_non_act,l_min,u_max;
    double norm_proj_d,sup_norm_proj_g,gd,q,q_old,gd_old,f_newton_first,delta_f0;
    std::vector<double> x,g,v,w,l,u,x_best,g_best,d,g_q,hp,lambda,mu,g1;
    bool hd_exact,act_phase,f_decreased,ls,ind_found,f_computed,gd_exit,dir_exit,stepsize_exit,min_f_exit;
    bool checkpoint,is_restarted,is_first_linesearch;
    bool warn_grad,warn_small,warn_noposdef,warn_sing,warn_conjfail,warn_maxcgit,is_d_neg_curv;
    std::ofstream file_output;

    bool converged(),check_stop(bool);
    void update_w(),restart();
    void trunc_newton_dir(),linesearch(),compute_sup_norm_proj_g(),approximate_hd(const std::vector<double>&);
    void estimate_active_set(),compute_multipliers();
    void main_prints();
    void err_obj(),clear_vectors();

public:
    Asa_bcp(unsigned short int&, Problem*, const asa_bcp_options*);
    Asa_bcp(unsigned short int& asa_bcp_status, Problem* p) : Asa_bcp(asa_bcp_status, p, NULL){};
    void solve();
    const std::vector<double>& get_x();
    const double get_f();
    const double get_sup_norm_proj_g();
    const unsigned int get_it();
    const unsigned int get_n_f();
    const unsigned int get_n_g();
    const unsigned int get_n_hd();
    const unsigned int get_inner_it();
    const short int get_flag();

};

#endif