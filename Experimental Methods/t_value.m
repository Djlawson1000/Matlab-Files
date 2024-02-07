function result = t_value(conf, nu)
    alpha = 1 - conf;
    pUp = 1 - alpha/2;
    result = tinv(pUp, nu);
end
