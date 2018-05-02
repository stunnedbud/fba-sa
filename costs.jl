
# Returns number of times solution v fails to meet given constraints
function count_constraints_fails(v::Array{Float64,2}, cons::Array{Float64,2})
    i = 1
    count = 0
    for x in v
        if x < cons[i,1] || x > cons[i,2]
            count += 1
        end
        i+=1
    end
    count
end

# Returns Σ_i=1 to len(z) of (z_i) , where z = S*v
# It quantifies how unstable a proposed solution is. A stable solution returns 0.
function stability_sum(v::Array{Float64,2}, S::Array{Float64,2})
    z = S * v
    sum(z)
end

# A solution v is feasible if it satisfies the equation S*v=0 and the constraints given
# for each element in v.
function feasible(v::Array{Float64}, S::Array{Float64,2}, cons::Array{Float64,2})
    # is the system stable? (within a certain range)
    if stability_sum(v,S)>1 || stability_sum(v,S) < -1 #stability_sum(v,S) != 0
        return false
    end
    
    # does it satisfy the given constraints? 
    if count_constraints_fails(v,cons) > 0
        return false
    end 
    
    true
end

# The cost function is simply biomass, which is given by the dot product v(dot)c
# There are two characteristics we wish to punish in a solution, the first is how many
# elements in v fail to meet constraints, and the second is how far away is the solution
# from being stable. This two values are pondered by the given punishment factor.
function cost(v::Array{Float64,2}, S::Array{Float64,2}, c::Array{Float64,2}, cons::Array{Float64,2}, pun::Float64)
    dot(v[:],c[:]) - pun * (count_constraints_fails(v,cons) + abs(stability_sum(v,S)))
end
