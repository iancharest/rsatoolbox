function f = portion(distribution_range, test_value)
    import rsa.util.*

    n = numel(distribution_range);
    n_lower = sumall(distribution_range < test_value);
    f = n_lower / n;
    
end%function
