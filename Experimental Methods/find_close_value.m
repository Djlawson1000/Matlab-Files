function index = find_close_value(value, array)
    array = double(array);  % Make sure this is a double
    [~, index] = min(abs(array - value));
end
