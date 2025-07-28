function posterior_appr = calculate_posterior(x_current,y_current,mu_x_appr,sd_x_appr,mu_y_appr,sd_y_appr,mu_x_av,sd_x_av,mu_y_av,sd_y_av,pi_appr,pi_av)

% Evaluate likelihoods
px_curr_appr = normpdf(x_current, mu_x_appr, sd_x_appr);
py_curr_appr = normpdf(y_current, mu_y_appr, sd_y_appr);
pXY_appr = px_curr_appr * py_curr_appr;  % naive independence

px_curr_av   = normpdf(x_current, mu_x_av,   sd_x_av);
py_curr_av   = normpdf(y_current, mu_y_av,   sd_y_av);
pXY_av   = px_curr_av * py_curr_av;

% Posterior probability state=1
posterior_appr = (pXY_appr * pi_appr) / ...
    (pXY_appr * pi_appr + pXY_av * pi_av);

end