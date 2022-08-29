#include <Rcpp.h>

#include "decision_tree.h"
#include "MCMC_utils.h"

using namespace Rcpp;

// [[Rcpp::export]]
List MCMC_sep(
    const NumericMatrix& Xpred,
    const NumericVector& Y_trt,
    NumericVector& Y_out,
    const double p_grow,   // Prob. of GROW
    const double p_prune,  // Prob. of PRUNE
    const double p_change, // Prob. of CHANGE
    const int m,           // Num. of Trees: default setting 100
    const int nu,
    double lambda, double lambda_1, double lambda_0,
    double dir_alpha, double alpha, double beta,
    const int n_iter,
    const bool verbose = false
) {

    // Data preparation
    const int P  = Xpred.ncol(); // number of covariates
    const int n  = Xpred.nrow(); // number of observations
    const int n1 = sum(Y_trt);   // number of treated sample
    const int n0 = n - n1;       // number of control sample

    Y_out = Y_out - mean(Y_out); // center mean

    NumericVector Y1(n1), Y0(n0);
    NumericMatrix Xpred1(n1, P), Xpred0(n0, P);
    {
        IntegerVector Y_trt_int = as<IntegerVector>(Y_trt);
        IntegerVector idx1(n1), idx0(n0);
        int count1 = 0, count0 = 0;
        for (int i = 0; i < n; i++)
        {
            if (Y_trt(i) == 1.0)
            {
                idx1(count1) = i; // idx of Y_trt = 1
                count1++;
            }
            else
            {
                idx0(count0) = i; // idx of Y_trt = 0
                count0++;
            }
        }
        for (int i = 0; i < n1; i++)
        {
            Xpred1(i, _) = Xpred(idx1(i), _); // potential confounders (A=1)
            Y1(i) = Y_out(idx1(i));           // Outcome (A=1)
        }
        for (int i = 0; i < n0; i++)
        {
            Xpred0(i, _) = Xpred(idx0(i), _); // potential confounders (A=0)
            Y0(i) = Y_out(idx0(i));           // Outcome (A=0)
        }
    }

    // Xcut <- lappy(1:dim(Xpred)[2], function(t) sort(unique(Xpred[,t])))
    NumericVector Xcut[P], Xcut1[P], Xcut0[P]; // e.g. unique value of potential confounders
    for (int j = 0; j < P; j++)
    {
        NumericVector temp;
        temp = unique(Xpred(_, j));
        temp.sort();
        Xcut[j] = clone(temp);

        temp = unique(Xpred1(_, j));
        temp.sort();
        Xcut1[j] = clone(temp);

        temp = unique(Xpred0(_, j));
        temp.sort();
        Xcut0[j] = clone(temp);
    }

    // Initial Setup
    // Priors, initial values and hyper-parameters
    NumericVector Z = Rcpp::rnorm(n, R::qnorm(mean(Y_trt), 0, 1, true, false), 1); // latent variable
    NumericVector prob = {p_grow, p_prune, p_change};

    double sigma2 = 1;
    NumericVector sigma2_1       (n_iter + 1); // create placeholder for sigma2_1
    NumericVector sigma2_0       (n_iter + 1); // create placeholder for sigma2_0
    NumericVector dir_alpha_hist (n_iter + 1); // create placeholder for dir_alpha
    sigma2_1(0)       = var(Y1);
    sigma2_0(0)       = var(Y0);
    dir_alpha_hist(0) = dir_alpha;

    // sigma_mu based on min/max of Y, Y (A=1) and Y (A=0)
    double sigma_mu   = std::max(pow(min(Z)  / (-2 * sqrt(m)), 2), pow(max(Z)  / (2 * sqrt(m)), 2));
    double sigma_mu_1 = std::max(pow(min(Y1) / (-2 * sqrt(m)), 2), pow(max(Y1) / (2 * sqrt(m)), 2));
    double sigma_mu_0 = std::max(pow(min(Y0) / (-2 * sqrt(m)), 2), pow(max(Y0) / (2 * sqrt(m)), 2));

    // Initial values of R
    NumericVector R  = clone(Z);
    NumericVector R1 = clone(Y1);
    NumericVector R0 = clone(Y0);

    // Initial values for the selection probabilities
    NumericVector post_dir_alpha  = rep(1.0, P);
    NumericVector post_dir_alpha1 = rep(1.0, P);
    NumericVector post_dir_alpha0 = rep(1.0, P);

    // thin = 10, burn-ins = n_iter/2
    int thin       = 10;
    int burn_in    = n_iter / 2;
    int n_post     = (n_iter - burn_in) / thin; // number of post sample
    int thin_count = 1;
    // IntegerVector seq = init_seq(n_iter, thin, burn_in);
    // IntegerVector::iterator seq_pt = seq.begin();

    NumericVector Effect (n_post);
    NumericVector PO_Y1  (n_post);
    NumericVector PO_Y0  (n_post);
    NumericMatrix predicted_Y1 (n, n_post);
    NumericMatrix predicted_Y0 (n, n_post);
    IntegerMatrix ind    (n_post, P);
    int post_sample_idx = 0;

    IntegerMatrix Obs_list  (n,  m); // changed list to matrix
    IntegerMatrix Obs1_list (n1, m);
    IntegerMatrix Obs0_list (n0, m);


    // Place-holder for the posterior samples
    NumericMatrix Tree   (n,     m);
    NumericMatrix Tree1  (n1,    m);
    NumericMatrix Tree0  (n0,    m);
    NumericMatrix Tree11 (n, m);
    NumericMatrix Tree00 (n, m);

    DecisionTree dt_list[m]; // changed to array of trees
    DecisionTree dt1_list[m];
    DecisionTree dt0_list[m];
    for (int t = 0; t < m; t++)
    {
        dt_list[t]  = DecisionTree(n,  t);
        dt1_list[t] = DecisionTree(n1, t);
        dt0_list[t] = DecisionTree(n0, t);
    }

    // Obtaining namespace of MCMCpack package
    Environment MCMCpack = Environment::namespace_env("MCMCpack");

    // Picking up rinvgamma() and rdirichlet() function from MCMCpack package
    Function rinvgamma  = MCMCpack["rinvgamma"];
    Function rdirichlet = MCMCpack["rdirichlet"];

    NumericVector prop_prob = rdirichlet(1, rep(dir_alpha, P));

    NumericMatrix xpred_mult(n, P + 1);                          // placeholder for bootstrap sample for inference
    IntegerVector xpred_ind = sample(Xpred.nrow(), n, true) - 1; // collecting indices of bootstrap samples
    for (int j = 0; j < P; j++)
    {
        for (int i = 0; i < n; i++)
        {
            xpred_mult(i, j) = Xpred(i, j);
        }
    }
    for (int i = 0; i < n; i++)
        xpred_mult(i, P) = i;

    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < P; j++)
            xpred_mult(i, j) = Xpred(i, j);
        xpred_mult(i, P) = i;
    }


    ////////////////////////////////////////
    //////////   Run main MCMC    //////////
    ////////////////////////////////////////

    // Run MCMC
    for (int iter = 1; iter <= n_iter; iter++)
    {
        if (verbose)
        {
            if (iter % 100 == 0)
                Rcout << "Rcpp iter : " << iter << " of " << n_iter << std::endl;
        }
        update_Z(Z, Y_trt, Tree);

        // ------ Exposure Model
        for (int t = 0; t < m; t++)
        {
            update_R(R, Z, Tree, t);

            if (dt_list[t].length() == 1)
            { 
                // tree has no node yet
                dt_list[t].GROW_first(
                    Xpred, Xcut, sigma2, sigma_mu, R, Obs_list,
                    p_prune, p_grow, alpha, beta, prop_prob);
            }
            else
            {
                int step = sample(3, 1, false, prob)(0);
                switch (step)
                {
                    case 1: // GROW step
                        dt_list[t].GROW(
                            Xpred, Xcut, sigma2, sigma_mu, R, Obs_list,
                            p_prune, p_grow, alpha, beta, prop_prob
                        );
                        break;

                    case 2: // PRUNE step
                        dt_list[t].PRUNE(
                            Xpred, Xcut, sigma2, sigma_mu, R, Obs_list,
                            p_prune, p_grow, alpha, beta, prop_prob
                        );
                        break;

                    case 3: // CHANGE step
                        dt_list[t].CHANGE(
                            Xpred, Xcut, sigma2, sigma_mu, R, Obs_list,
                            p_prune, p_grow, alpha, beta, prop_prob
                        );
                        break;

                    default: {};
                } // end of switch
            }     // end of tree instance
            dt_list[t].Mean_Parameter(Tree, sigma2, sigma_mu, R, Obs_list);
        } // end of Exposure Model

        // ------ Outcome Model (A=1)
        for (int t = 0; t < m; t++)
        {
            update_R(R1, Y1, Tree1, t);
            if (dt1_list[t].length() == 1)
            { 
                // tree has no node yet
                dt1_list[t].GROW_first(
                    Xpred1, Xcut1, sigma2_1(iter - 1), sigma_mu_1, R1, Obs1_list,
                    p_prune, p_grow, alpha, beta, prop_prob
                );
            }
            else
            {
                int step = sample(3, 1, false, prob)(0);
                switch (step)
                {
                    case 1: // GROW step
                        dt1_list[t].GROW(
                            Xpred1, Xcut1, sigma2_1(iter - 1), sigma_mu_1, R1, Obs1_list,
                            p_prune, p_grow, alpha, beta, prop_prob
                        );
                        break;

                    case 2: // PRUNE step
                        dt1_list[t].PRUNE(
                            Xpred1, Xcut1, sigma2_1(iter - 1), sigma_mu_1, R1, Obs1_list,
                            p_prune, p_grow, alpha, beta, prop_prob
                        );
                        break;

                    case 3: // CHANGE step
                        dt1_list[t].CHANGE(
                            Xpred1, Xcut1, sigma2_1(iter - 1), sigma_mu_1, R1, Obs1_list,
                            p_prune, p_grow, alpha, beta, prop_prob
                        );
                        break;

                    default: {};
                } // end of switch
            }     // end of tree instance
            dt1_list[t].Mean_Parameter(Tree1, sigma2_1(iter - 1), sigma_mu_1, R1, Obs1_list);
        } // end of Outcome Model (A=1)

        // ------ Outcome Model (A=0)
        for (int t = 0; t < m; t++)
        {
            update_R(R0, Y0, Tree0, t);
            if (dt0_list[t].length() == 1)
            { 
                // tree has no node yet
                dt0_list[t].GROW_first(
                    Xpred0, Xcut0, sigma2_0(iter - 1), sigma_mu_0, R0, Obs0_list,
                    p_prune, p_grow, alpha, beta, prop_prob
                );
            }
            else
            {
                int step = sample(3, 1, false, prob)(0);
                switch (step)
                {
                    case 1: // GROW step
                        dt0_list[t].GROW(
                            Xpred0, Xcut0, sigma2_0(iter - 1), sigma_mu_0, R0, Obs0_list,
                            p_prune, p_grow, alpha, beta, prop_prob
                        );
                        break;

                    case 2: // PRUNE step
                        dt0_list[t].PRUNE(
                            Xpred0, Xcut0, sigma2_0(iter - 1), sigma_mu_0, R0, Obs0_list,
                            p_prune, p_grow, alpha, beta, prop_prob
                        );
                        break;

                    case 3: // CHANGE step
                        dt0_list[t].CHANGE(
                            Xpred0, Xcut0, sigma2_0(iter - 1), sigma_mu_0, R0, Obs0_list,
                            p_prune, p_grow, alpha, beta, prop_prob
                        );
                        break;

                    default: {};
                } // end of switch
            }     // end of tree instance
            dt0_list[t].Mean_Parameter(Tree0, sigma2_0(iter - 1), sigma_mu_0, R0, Obs0_list);
        } // end of Outcome Model (A=0)

        //  Sample variance parameter
        {
            NumericVector sigma2_1_temp = rinvgamma(1, nu / 2 + n1 / 2, nu * lambda_1 / 2 + sum(pow(Y1 - rowSums(Tree1), 2)) / 2);
            NumericVector sigma2_0_temp = rinvgamma(1, nu / 2 + n0 / 2, nu * lambda_0 / 2 + sum(pow(Y0 - rowSums(Tree0), 2)) / 2);

            sigma2_1(iter) = sigma2_1_temp(0);
            sigma2_0(iter) = sigma2_0_temp(0);
        }

        // Num. of inclusion of each potential confounder
        NumericVector add(P), add1(P), add0(P);
        for (int t = 0; t < m; t++)
        {
            add  += dt_list[t].num_included(P);
            add1 += dt1_list[t].num_included(P);
            add0 += dt0_list[t].num_included(P);
        }

        // M.H. algorithm for the alpha parameter in the dirichlet distribution (after some warm-ups)
        if (iter < n_iter / 10)
        {
            post_dir_alpha = rep(1.0, P) + add + add1 + add0;
        }
        else
        {
            double p_dir_alpha = std::max(R::rnorm(dir_alpha, 0.1), pow(0.1, 10));

            NumericVector SumS(P);
            log_with_LB(SumS, prop_prob);

            double dir_lik_p, dir_lik, ratio;

            dir_lik_p =
                sum(SumS * (rep(p_dir_alpha / P, P) - 1)) + lgamma(sum(rep(p_dir_alpha / P, P))) - sum(lgamma(rep(p_dir_alpha / P, P)));

            dir_lik =
                sum(SumS * (rep(dir_alpha / P, P) - 1)) + lgamma(sum(rep(dir_alpha / P, P))) - sum(lgamma(rep(dir_alpha / P, P)));

            ratio =
                dir_lik_p + log(pow(p_dir_alpha / (p_dir_alpha + P), 0.5 - 1) * pow(P / (p_dir_alpha + P), 1 - 1) * abs(1 / (p_dir_alpha + P) - p_dir_alpha / pow(p_dir_alpha + P, 2))) + R::dnorm(dir_alpha, p_dir_alpha, 0.1, true) - dir_lik - log(pow(dir_alpha / (dir_alpha + P), 0.5 - 1) * pow(P / (dir_alpha + P), 1 - 1) * abs(1 / (dir_alpha + P) - dir_alpha / pow(dir_alpha + P, 2))) - R::dnorm(p_dir_alpha, dir_alpha, 0.1, true);

            if (ratio > log(R::runif(0, 1)))
            {
                dir_alpha = p_dir_alpha;
            }

            post_dir_alpha = rep(dir_alpha / P, P) + add + add1 + add0;
        } // end of M.H. algorithm

        dir_alpha_hist(iter) = dir_alpha;
        prop_prob = rdirichlet(1, post_dir_alpha);

        // Sampling E[Y(1)-Y(0)]
        if (iter > burn_in)
        {
            if (thin_count < thin)
            {
                thin_count++;
            }
            else
            {
                thin_count = 1;
                for (int i = 0; i < m; i++)
                {
                    dt1_list[i].Predict_sep(Tree11, Xcut1, xpred_mult, n);
                    dt0_list[i].Predict_sep(Tree00, Xcut0, xpred_mult, n);
                }
                // Effect(post_sample_idx) = mean(rowSums(Tree11) - rowSums(Tree00));
                predicted_Y1 (_, post_sample_idx) = clone(rowSums(Tree11));
                predicted_Y0 (_, post_sample_idx) = clone(rowSums(Tree00));
                PO_Y1  (post_sample_idx)          = mean(predicted_Y1(_, post_sample_idx));
                PO_Y0  (post_sample_idx)          = mean(predicted_Y0(_, post_sample_idx));
                Effect (post_sample_idx)          = PO_Y1(post_sample_idx) - PO_Y0(post_sample_idx);

                IntegerVector ind_temp  = ifelse(add1 + add0 > 0.0, 1, 0);
                ind(post_sample_idx, _) = clone(ind_temp); // indicator of whether confounders are included
                post_sample_idx++;
            }
        }

        Rcpp::checkUserInterrupt(); // check for break in R
    } // end of MCMC iterations

    List L = List::create(
        Named("Effect")       = Effect,
        Named("PO_Y1")        = PO_Y1,
        Named("PO_Y0")        = PO_Y0,
        Named("predicted_Y1") = predicted_Y1,
        Named("predicted_Y0") = predicted_Y0,
        Named("xpred_mult")   = xpred_mult,
        Named("ind")          = ind,
        Named("sigma2_1")     = sigma2_1,
        Named("sigma2_0")     = sigma2_0,
        Named("dir_alpha")    = dir_alpha_hist
    );

    return L;
}
