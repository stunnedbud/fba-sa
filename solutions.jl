# Generates a random solution within the given constraints for elements of v.
function random_solution(cons::Array{Float64,2})
    v = fill(0.0, (length(cons)/2, 1))
    for i in 1:length(cons)
        v[i] = rand(cons[i,1]:cons[i,2])
    end
    v
end

# Adds/subtracts 1 to a random element of v. Returns new array.
function neighbor(v::Array{Float64,2})
    i = rand(1:length(v))
    u = copy(v)
    if rand() > 0.5 
        u[i] += 1
    else
        u[i] -= 1
    end
    u
end
