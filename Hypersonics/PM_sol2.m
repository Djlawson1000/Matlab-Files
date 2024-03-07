function [M2,P2,T2] = PM_sol2(theta1,M1,P1,T1,gamma)
    
    theta1 = deg2rad(theta1); 

    gg = (gamma + 1) / (gamma - 1);
    vPM = @(Md) (sqrt(gg) * atan(sqrt(1/gg*(Md^2-1))) - atan(sqrt(Md^2-1))) -...
        (sqrt(gg) * atan(sqrt(1/gg*(M1^2-1))) - atan(sqrt(M1^2-1))) - theta1;
    M2 = fzero(vPM,M1,optimset('Display','off'));
    
    gg = (gamma - 1)/2;
    ggg = gamma/(gamma - 1);
    T2 = T1 * (1 + gg * M1^2) / (1 + gg * M2^2);
    P2 = P1 * (((1 + gg * M1^2) / (1 + gg * M2^2))^ggg);

end