
# Generates a random solution within the given constraints for elements of v.
function random_solution(cons::Array{Float64,2})
    vlen = convert(Int, floor(length(cons)/2))
    v = fill(0.0, (vlen, 1))
    for i in 1:vlen
        v[i] = rand(cons[i,1]:cons[i,2])
    end
    v
end

# 
function linear_combination(B::Array{Float64,2}, c::Array{Float64,1})
    v = fill(0.0, (length(B[:,1]),1))
    for i in 1:length(c)
        v += c[i]*B[:,i]
    end
    v
end

# Returns basis for kernel of S, as well as a random vector formed by a linear
# combination of said basis.
function stable_solution(B::Array{Float64,2}, region::Float64, range::Float64)
    # get basis for nullspace of S
    # each of the columns of B is a vector b that satisfies Sb=0
    # also any linear combination of them satisfies it
    #B = nullspace(S) ## turns out its cheaper to compute once and pass it everywhere
    
    # Generate random vector of scalars to combine the basis vectors with
    vectors_in_basis = length(B[1,:])
    c = fill(region, vectors_in_basis)
    for i in 1:vectors_in_basis
        a = rand(0:range)
        if rand() > 0.5
            c[i] += a
        else
            c[i] -= a
        end
    end
    
    # Combine vectors in basis to form solution v
    v = linear_combination(B,c)
    c, v
end

# Returns array with indexes of elements in v that fail constraints cons
function failed_constraints(v::Array{Float64,2}, cons::Array{Float64,2})
    r = Int[]
    for i in 1:length(v)
        if v[i] < cons[i,1] || v[i] > cons[i,2]
            append!(r,i)
        end
    end
    r
end

#
function random_failed_constraint(v::Array{Float64,2}, cons::Array{Float64,2})
    return -1
end

#
function neighbor(c::Array{Float64,1}, v::Array{Float64,2}, cons::Array{Float64,2}, B::Array{Float64,2}, step_size::Float64)
    R = failed_constraints(v,cons)
    i = rand(1:length(c))
    
    d = copy(c) # neighbor
    v2 = Union{}
   
    # tests
    #x = step_size * B[:,i]
    #v3 = v + x
    #d[i] += step_size
    #v4 = linear_combination(B,d)
    #for j in 1:95
    #    println(v[j])
    #    println(x[j])
    #    println(v3[j])    
    #    println(v4[j])
    #    println(" ")
    #end
    #println("DUN")
    
    
    b = step_size * B[:,i] # random vector in basis times step_size
    if rand() > 0.5
        v2 = v + b
        d[i] += step_size
    else
        v2 = v - b
        d[i] -= step_size
    end
    v2 = linear_combination(B,d)
    d,v2
end

# Adds/subtracts 1 to a random element of v, staying within the given constraints. 
# Returns new array.
#=
function neighbor(v::Array{Float64,2}, cons::Array{Float64,2}, step_size::Float64)
    i = rand(1:length(v))
    u = copy(v)
    if rand() > 0.5 
        if u[i] + 1 > cons[i,2]
            u[i] -= step_size
        else
            u[i] += step_size
        end
    else
        if u[i] -1 < cons[i,1]
            u[i] += step_size
        else
            u[i] -= step_size
        end
    end
    u
end
=#
