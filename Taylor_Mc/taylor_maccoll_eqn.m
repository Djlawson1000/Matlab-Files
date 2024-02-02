%this function solves the state vairble form of the Taylor-Maccoll Eq. 

function [value, isterminal, direction] = taylor_maccoll_eqn(t,g,flag)
G = 1.4;

gdot(1)=g(2);

a = 0.5*(G - 1) ;

Numerator = (a.*(2*g(1) + g(2).*cot(t) - 2*(g(1).^3) - (g(1).^2).*(g(2)).*cot(t) - (2*(g(1)).*(g(2).^2)) - (cot(t).*g(2).^3))) - (g(1).*(g(2).^2));
Denominator =  g(2).^2 + a.*((g(1).^2) + (g(2).^2) - 1);
gdot(2) = Numerator./Denominator;
% this event detects when g(2)(angular velocity) is zero, stop the iteration from proceeding

if (nargin<3)||isempty(flag)

value=[gdot(1);gdot(2)];

%used youtube video(reference 3) to understand the functioning of ode23
else flag=='events';
value=g;            % values to check
isterminal=[0;1];   % terminal on g(2)
direction=[0;1];   % detects falling slope(not sure what this does)
end

%end of subroutine