function [sig_fit, eps_fit] = polyfit_err(x, y, order, p, conf)
    N = length(x);
    nu = N - (order + 1);
    alpha = 1 - conf;
    pUp = 1 - alpha / 2;
    t_val = tinv(pUp, nu);
    
    y_new = polyval(p, x);
    
    sig_fit = sqrt(sum((y - y_new).^2) / nu);
    
    eps_fit = t_val * sig_fit;
end
