function result = round_to_nsf(number, nsf)
    if nargin < 2
        nsf = 2;
    end
    
    integer_part = floor(number);
    result = round(number, nsf - numel(num2str(integer_part)));
end
