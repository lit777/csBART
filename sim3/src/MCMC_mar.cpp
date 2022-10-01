#include <Rcpp.h>

#include "decision_tree.h"
#include "MCMC_utils.h"

using namespace Rcpp;

// [[Rcpp::export]]
List MCMC_mar(
    const NumericMatrix& Xpred,
    const NumericVector& Y_trt,
    NumericVector& Y_out,
    const double p_grow,   // Prob. of GROW
    const double p_prune,  // Prob. of PRUNE
    const double p_change, // Prob. of CHANGE
    const int m,           // Num. of Trees: default setting 100
    const int nu,
    double lambda, double lambda_1,
    double dir_alpha, double alpha, double beta,
    const int n_iter,
    const bool verbose = false
) {

    // Data preparation
    const int P = Xpred.ncol(); // number of covariates
    const int n = Xpred.nrow(); // number of observations

    NumericMatrix Xpred1 = cbind(Y_trt, Xpred); 

    NumericVector Xcut[P], Xcut1[P + 1]; // e.g. unique value of potential confounders
    for (int j = 0; j < P; j++)
    {
        NumericVector temp = unique(Xpred(_, j));
        temp.sort();
        Xcut[j] = temp;
    }
    for (int j = 0; j < P + 1; j++)
    {
        NumericVector temp = unique(Xpred1(_, j));
        temp.sort();
        Xcut1[j] = temp;
    }


    // Initial Setup
    // Priors, initial values and hyper-parameters
    NumericVector Z = Rcpp::rnorm(n, R::qnorm(mean(Y_trt), 0, 1, true, false), 1); // latent variable
    NumericVector prob = {p_grow, p_prune, p_change};

    double sigma2 = 1;
    NumericVector sigma2_1 (n_iter + 1);       // create placeholder for sigma2_1
    NumericVector dir_alpha_hist (n_iter + 1); // create placeholder for dir_alpha
    sigma2_1(0)       = var(Y_out);
    dir_alpha_hist(0) = dir_alpha;

    // sigma_mu based on min/max of Y, Y (A=1) and Y (A=0)
    double sigma_mu   = std::max(pow(min(Z)     / (-2 * sqrt(m)), 2), pow(max(Z)     / (2 * sqrt(m)), 2));
    double sigma_mu_1 = std::max(pow(min(Y_out) / (-2 * sqrt(m)), 2), pow(max(Y_out) / (2 * sqrt(m)), 2));

    // Initial values of R
    NumericVector R  = clone(Z);
    NumericVector R1 = clone(Y_out);

    // Initial values for the selection probabilities
    NumericVector post_dir_alpha  = rep(1.0, P);
    NumericVector post_dir_alpha1 = rep(1.0, P + 1);

    // thin = 10, burn-ins = n_iter/2
    int thin       = 10;
    int burn_in    = n_iter / 2;
    int n_post     = (n_iter - burn_in) / thin; // number of post sample
    int thin_count = 1;

    NumericVector Effect (n_post);
    NumericVector PO_Y1  (n_post);
    NumericVector PO_Y0  (n_post);
    NumericMatrix predicted_Y1 (n, n_post);
    NumericMatrix predicted_Y0 (n, n_post);
    IntegerMatrix ind    (n_post, P+1);
    int post_sample_idx = 0;

    IntegerMatrix Obs_list(n, m); // changed list to matrix
    IntegerMatrix Obs1_list(n, m);


    // Place-holder for the posterior samples
    NumericMatrix Tree   (n,     m);
    NumericMatrix Tree1  (n,     m);
    NumericMatrix Tree11 (n, m);
    NumericMatrix Tree00 (n, m);

    DecisionTree dt_list[m]; // changed to array of trees
    DecisionTree dt1_list[m];
    for (int t = 0; t < m; t++)
    {
        dt_list[t]  = DecisionTree(n, t, true, false);
        dt1_list[t] = DecisionTree(n, t, false, true);
    }

    // Obtaining namespace of MCMCpack package
    Environment MCMCpack = Environment::namespace_env("MCMCpack");

    // Picking up rinvgamma() and rdirichlet() function from MCMCpack package
    Function rinvgamma  = MCMCpack["rinvgamma"];
    Function rdirichlet = MCMCpack["rdirichlet"];

    NumericVector prop_prob = rdirichlet(1, rep(dir_alpha, P + 1));

    NumericMatrix xpred_mult(n, P + 2);                          // placeholder for bootstrap sample for inference
    IntegerVector xpred_ind = sample(Xpred.nrow(), n, true) - 1; // collecting indices of bootstrap samples
    for (int j = 1; j < P + 1; j++)
    { 
        // first col is rep(0, 2n)
        for (int i = 0; i < n; i++)
        {
            xpred_mult(i, j) = Xpred1(i, j);
        }
    }
    for (int i = 0; i < n; i++)
        xpred_mult(i, P + 1) = i;

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
            // decision trees
            update_R(R, Z, Tree, t);

            NumericVector prop_prob_exp = prop_prob[Rcpp::Range(1, P)];
            prop_prob_exp = prop_prob_exp / (sum(prop_prob) - prop_prob(0));

            if (dt_list[t].length() == 1)
            { 
                // tree has no node yet
                // grow first step
                dt_list[t].GROW_first(
                    Xpred, Xcut, sigma2, sigma_mu, R, Obs_list,
                    p_prune, p_grow, alpha, beta, prop_prob_exp
                );
            }
            else
            {
                int step = sample(3, 1, false, prob)(0);
                switch (step)
                {
                    case 1: // GROW step
                        dt_list[t].GROW(
                            Xpred, Xcut, sigma2, sigma_mu, R, Obs_list,
                            p_prune, p_grow, alpha, beta, prop_prob_exp
                        );
                        break;

                    case 2: // PRUNE step
                        dt_list[t].PRUNE(
                            Xpred, Xcut, sigma2, sigma_mu, R, Obs_list,
                            p_prune, p_grow, alpha, beta, prop_prob_exp
                        );
                        break;

                    case 3: // CHANGE step
                        dt_list[t].CHANGE(
                            Xpred, Xcut, sigma2, sigma_mu, R, Obs_list,
                            p_prune, p_grow, alpha, beta, prop_prob_exp
                        );
                        break;

                    default: {}
                } // end of switch
            }     // end of tree instance
            dt_list[t].Mean_Parameter(Tree, sigma2, sigma_mu, R, Obs_list);
        } // end of Exposure Model

        // ------ Outcome Model
        for (int t = 0; t < m; t++)
        {
            update_R(R1, Y_out, Tree1, t);
            if (dt1_list[t].length() == 1)
            { 
                // tree has no node yet
                // grow_first step
                dt1_list[t].GROW_first(
                    Xpred1, Xcut1, sigma2_1(iter - 1), sigma_mu_1, R1, Obs1_list,
                    p_prune, p_grow, alpha, beta, prop_prob);
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

                    default: {}
                } // end of switch
            }     // end of tree instance
            dt1_list[t].Mean_Parameter(Tree1, sigma2_1(iter - 1), sigma_mu_1, R1, Obs1_list);
        } // end of Outcome Model (A=1)

        // Sample variance parameter
        {
            // sigma2 = 1;
            NumericVector sigma2_1_temp = rinvgamma(1, nu / 2 + n / 2, nu * lambda_1 / 2 + sum(pow(Y_out - rowSums(Tree1), 2)) / 2);
            sigma2_1(iter) = sigma2_1_temp(0);
        }

        // Num. of inclusion of each potential confounder
        NumericVector add(P), add1(P + 1);
        for (int t = 0; t < m; t++)
        {
            add  += dt_list[t].num_included(P);
            add1 += dt1_list[t].num_included(P + 1);
        }

        // M.H. algorithm for the alpha parameter in the dirichlet distribution
        {
            double p_dir_alpha, dir_lik_p, dir_lik, ratio;

            p_dir_alpha = std::max(R::rnorm(dir_alpha, 0.1), pow(0.1, 10));

            NumericVector SumS(P + 1);
            log_with_LB(SumS, prop_prob);

            dir_lik_p =
                sum(SumS * (rep(p_dir_alpha / (P + 1), (P + 1)) - 1)) + lgamma(sum(rep(p_dir_alpha / (P + 1), (P + 1)))) - sum(lgamma(rep(p_dir_alpha / (P + 1), (P + 1))));

            dir_lik =
                sum(SumS * (rep(dir_alpha / (P + 1), (P + 1)) - 1)) + lgamma(sum(rep(dir_alpha / (P + 1), (P + 1)))) - sum(lgamma(rep(dir_alpha / (P + 1), (P + 1))));

            ratio = dir_lik_p + log(pow(p_dir_alpha / (p_dir_alpha + (P + 1)), 0.5 - 1) * abs(1 / (p_dir_alpha + (P + 1)) - p_dir_alpha / pow(p_dir_alpha + (P + 1), 2))) + R::dnorm(dir_alpha, p_dir_alpha, 0.1, true) - dir_lik - log(pow(dir_alpha / (dir_alpha + (P + 1)), 0.5 - 1) * abs(1 / (dir_alpha + (P + 1)) - dir_alpha / pow(dir_alpha + (P + 1), 2))) - R::dnorm(p_dir_alpha, dir_alpha, 0.1, true);

            if (ratio > log(R::runif(0, 1)))
            {
                dir_alpha = p_dir_alpha;
            }
            dir_alpha_hist(iter) = dir_alpha;
            post_dir_alpha = rep(dir_alpha / (P + 1), P + 1);
        }

        // M.H. algorithm for the inclusion probabilities
        {
            double dir_lik_p, dir_lik, ratio;

            NumericVector add_max(add.length() + 1);
            add_max(0) = add1(0);
            add_max[Rcpp::Range(1, add.length())] = add;

            NumericVector p_prop_prob = rdirichlet(1, add1 + post_dir_alpha + add_max);
            NumericVector log_p_prop_prob(P + 1), log_prop_prob(P + 1);
            log_with_LB(log_p_prop_prob, p_prop_prob);
            log_with_LB(log_prop_prob,   prop_prob);

            dir_lik_p = sum(add) * log(1 / (1 - p_prop_prob(0))) + (add1(0) + post_dir_alpha(0) - 1.0) * log_p_prop_prob(0) + sum((add1[Rcpp::Range(1, P)] + add + post_dir_alpha[Rcpp::Range(1, P)] - 1.0) * log_p_prop_prob[Rcpp::Range(1, P)]);

            dir_lik = sum(add) * log(1 / (1 - prop_prob(0))) + (add1(0) + post_dir_alpha(0) - 1.0) * log_prop_prob(0) + sum((add1[Rcpp::Range(1, P)] + add + post_dir_alpha[Rcpp::Range(1, P)] - 1.0) * log_prop_prob[Rcpp::Range(1, P)]);

            ratio = dir_lik_p + sum((add1 + post_dir_alpha + add_max - 1.0) * log_prop_prob) - dir_lik - sum((add1 + post_dir_alpha + add_max - 1.0) * log_p_prop_prob);

            if (ratio > log(R::runif(0, 1)))
            {
                prop_prob = clone(p_prop_prob);
            }
        }

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
                    dt1_list[i].Predict_mar(Tree11, Xcut1, xpred_mult, n, 1);
                    dt1_list[i].Predict_mar(Tree00, Xcut1, xpred_mult, n, 0);
                }
                // Effect(post_sample_idx) = mean(rowSums(Tree11) - rowSums(Tree00));
                predicted_Y1 (_, post_sample_idx) = clone(rowSums(Tree11));
                predicted_Y0 (_, post_sample_idx) = clone(rowSums(Tree00));
                PO_Y1  (post_sample_idx)          = mean(predicted_Y1(_, post_sample_idx));
                PO_Y0  (post_sample_idx)          = mean(predicted_Y0(_, post_sample_idx));
                Effect (post_sample_idx)          = PO_Y1(post_sample_idx) - PO_Y0(post_sample_idx);


                IntegerVector ind_temp  = ifelse(add1 > 0.0, 1, 0);
                ind(post_sample_idx, _) = ind_temp; // indicator of whether confounders are included
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
        Named("dir_alpha")    = dir_alpha_hist
    );

    return L;
}
