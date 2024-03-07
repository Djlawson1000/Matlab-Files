function [M2,P2,T2] = oblique2(theta1,M1,P1,T1,gamma)

    beta_out = wave_angle2(M1,theta1,gamma);
    beta_out = deg2rad(beta_out);

    theta1 = deg2rad(theta1);
    M1_norm = M1*sin(beta_out);

    M2_norm = (1+0.5*(gamma-1)*M1_norm^2)/(gamma*M1_norm^2-0.5*(gamma-1));
    M2_norm = sqrt(M2_norm);
    M2 = M2_norm/sin(beta_out-theta1);
    
    P2 = P1*((1+gamma*M1_norm^2)/(1+gamma*M2_norm^2));
    T2 = T1*((1+0.5*(gamma-1)*M1_norm^2)/(1+0.5*(gamma-1)*M2_norm^2));

end

