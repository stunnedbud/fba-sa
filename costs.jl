
# A solution v is feasible if it satisfies the equation S*v=0 and the constraints given
# for each element in v.
function feasible(v::Array{Float64}, S::Array{Float64,2}, cons::Array{Float64,2})
    # is the system stable?
    z = S * v
    for i in z
        if i != 0
            return false
        end
    end
    
    # does it satisfy the given constraints? 
    if count_constraints_fails(v,cons) > 0
        return false
    end 
    
    true
end

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

# The cost function is simply biomass, which is given by the dot product v(dot)c
function cost(v::Array{Float64,2}, c::Array{Float64,2}, cons::Array{Float64,2}, pun::Float64)
    dot(v[:],c[:]) - pun * count_constraints_fails(v,cons)
end
