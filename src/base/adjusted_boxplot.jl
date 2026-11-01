# Conclusion for outlier detection among marker-specific error rates : 
#  logit/log transformations and adjusted boxplot resulted in too few outliere
#  3IQR rule is the best choice since there are ususally many markers and a few false detections are thus acceptable


# If the adjusted boxplot is hiding points you believe are truly outliers, consider these "middle ground" approaches in Julia:
# The 3.0x IQR Rule: Instead of the standard 1.5, use the "Tukey Far-Out" multiplier of 3.0. 
# It’s a constant, but it’s often more intuitive than the exponential MC scaling.


# doi:10.1016/j.csda.2007.11.008
# An adjusted boxplot for skewed distributions
# using Statistics

function medcouple(x::AbstractVector)
    x = sort(x)
    # n = length(x)
    m = median(x)
    
    # Scale data around median
    z = x .- m
    z⁺ = z[z .>= 0]
    z⁻ = z[z .<= 0]
    p = length(z⁺)
    q = length(z⁻)
    
    # Kernel function h(z⁺_i, z⁻_j)
    # Handles special case where both are 0
    h(i, j) = if z⁺[i] == z⁻[j]
        sign(p - 1 - i - j)
    else
        (z⁺[i] + z⁻[j]) / (z⁺[i] - z⁻[j])
    end
    
    # Calculate all h values for the median (O(n²) implementation)
    h_values = [h(i, j) for i in 1:p, j in 1:q]
    return median(h_values)
end

function adjusted_fences(data)
    q1, q3 = quantile(data, [0.25, 0.75])
    iqr = q3 - q1
    mc = medcouple(data)
    
    if mc >= 0
        lower = q1 - 1.5 * exp(-4 * mc) * iqr
        upper = q3 + 1.5 * exp(3 * mc) * iqr
    else
        lower = q1 - 1.5 * exp(-3 * mc) * iqr
        upper = q3 + 1.5 * exp(4 * mc) * iqr
    end
    
    return (lower, upper, mc)
end
