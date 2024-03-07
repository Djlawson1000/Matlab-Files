function [beta_out] = wave_angle2(M1,theta1,gamma)
    
    theta1 = deg2rad(theta1);

    beta_guess = asin(1/M1);
    Func = @(beta) tan(theta1) - ...
        ((M1^2*sin(2*beta) - 2*cot(beta))/(M1^2*(gamma+cos(2*beta)) + 2));
    beta_out = fsolve(Func,beta_guess,optimset('Display', 'off'));
    
    beta_out = rad2deg(beta_out);
end

