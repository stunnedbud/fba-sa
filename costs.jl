
# Returns number of times solution v fails to meet given constraints
function count_constraints_fails(v::Array{Float64,2}, cons::Array{Float64,2}, precision::Float64)
    f = constraints_fails(v,cons,precision)
    length(f)
end

# Returns list with indices of elements of v that fail the given constraints
# (within the margin of error specified by the parameter precision)
function constraints_fails(v::Array{Float64,2}, cons::Array{Float64,2}, precision::Float64)
    i = 1
    fails = Int[]
    for x in v
        if x + precision < cons[i,1] || x - precision > cons[i,2]
            append!(fails, i) 
        end
        i+=1
    end
    fails
end

# Returns [Σ_i=1 to len(z) of (|z_i|)], where z = S*v
# It quantifies how unstable a proposed solution is. A perfectly stable solution returns 0.
# An acceptable solution (due to numerical stability issues) is some decimals off (around 10^-14)
function stability_sum(v::Array{Float64,2}, S::Array{Float64,2})
    z = S * v
    sum(abs.(z))
end

# A solution v is feasible if it satisfies the equation S*v=0 and the constraints given
# for each element in v.
function feasible(v::Array{Float64}, S::Array{Float64,2}, cons::Array{Float64,2}, precision::Float64)
    # is the system stable? (within a certain range)
    if stability_sum(v,S)>precision #stability_sum(v,S) != 0
        return false
    end
    
    # does it satisfy the given constraints? 
    if count_constraints_fails(v,cons,precision) > 0
        return false
    end 
    
    true
end

# There are two characteristics we wish to punish in a solution, the first is how many
# elements in v fail to meet constraints, and the second is how far away is the solution
# from being stable. This two values are pondered by the given punishment factor.
# If neither punishments apply, the cost function is simply to evaluate the objective 
# function, which is given by the dot product v(dot)c.
function cost(v::Array{Float64,2}, S::Array{Float64,2}, c::Array{Float64,2}, cons::Array{Float64,2}, pun1::Float64, pun2::Float64, pun3::Float64, precision::Float64)#, pun2::Float64)
    #pun2 = 200
    #pun3 = 1000
    #dot(v[:],c[:]) - pun * (count_constraints_fails(v,cons) + stability_sum(v,S))
    # since the solutions we are generating already satisfy Sv=0 there's no need to consider it here 
    
    #dot(v[:],c[:]) - pun * count_constraints_fails(v,cons) 
    
    # Quantify how far the elements that fail constraints are from satisfying them
    s = 0
    f = constraints_fails(v,cons,precision)
    for i in f
        if v[i,1] < cons[i,1] 
            s += cons[i,1] - v[i,1] #+ precision
        else # v[i] > cons[i,2]
            s += v[i,1] - cons[i,2] #- precision
        end
    end
    
    if length(f) <= 0 # we satisfy almost all constraints, start trying to maximize objective function
        return pun3 + (dot(v[:],c[:])) #- pun1 * s
        #if length(f) == 0 
        #    return pun3 + (dot(v[:],c[:])) #- pun1 * s 
        #else
        #    return 0 - pun1 * s  
        #end
    end # else ...
    0 - ( pun1 * s  + pun2 * length(f))
end
