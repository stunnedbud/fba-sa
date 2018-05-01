include("metabolic_network.jl")

function calc_lot(m::Array{Float64,2}, v::Array{Float64,2}, c::Array{Float64,2}, cons::Array{Float64,2}, T::Float64, L::Int, pun::Float64, plot_name::Int)
    v1 = v
    best = v
    cost_best = cost(v, c, cons, pun)
    accepted_costs = [cost_best]
    count = 0
    attempts = 0
    max_attempts = L^2
    cost_v1 = cost(v1, c, cons, pun)
    while count < L
        attempts += 1
        cost_v1 = cost(v1, c, cons, pun)
        v2 = neighbor(v1)
        cost_v2 = cost(v2, c, cons, pun)
        if cost_v2 >= cost_v1 - T # accepts solution
            v1 = v2
            cost_v1 = cost_v2
            count += 1
            push!(accepted_costs, cost_v1)
        end
        if cost_v2 > cost_best
            best = v2
            cost_best = cost_v2
        end
        if attempts > max_attempts
            println("Exceeded max attempts")
            return sum(accepted_costs)/L, v1, best, true # true to end run
        end 
    end
    
    # Uncomment to save timeseries of accepted costs to file.
    # Deactivated it because my memory was running out.
    out = ""
    for c in accepted_costs
        out *= string(c)*"\n" 
    end
    open("results/plots/"*string(plot_name)*".txt","a") do f
        write(f, out)
    end

    sum(accepted_costs)/L, v1, best, false
end


# "Main" function in this file. Calculates lot until average of accepted solutions improves, whereupon it
# decreases T by multiplying it by theta (cooling factor). Keeps doing that until T reaches epsilon.
# Returns last accepted solution and the absolute best found.
function acceptance_by_thresholds(m::Array{Float64,2}, v::Array{Float64,2}, c::Array{Float64,2}, cons::Array{Float64,2}, T::Float64, L::Int, ε::Float64, ϕ::Float64, pun::Float64, plot_name::Int)
    open("results/plots/"*string(plot_name)*".txt","w") do f
        write(f, "")
    end
    #srand(seed)
    p = 0
    abs_b = v # absolute best
    b = v
    f = false # flag for exceeding max attempts
    while T > ε
        q = p - 1
        while p > q
            q = p
            p, v, b, f = calc_lot(m, v, c, cons, T, L, pun, plot_name)    

            if cost(b,c,cons,pun) > cost(abs_b,c,cons,pun)
                abs_b = b
            end
            
            if f # exceeded max attempts
                println("FUCK")
                #b = sweep(g,b,pun,avg)
                return v,abs_b
            end
            
            # Sweeping
            #if cost(g,b,pun,avg) < cost(g,abs_b,pun,avg)
                #println("SWEEP:")
                #sw = sweep(g,b,pun,avg)
                #if cost(g,sw,pun,avg) < cost(g,abs_b,pun,avg)
                #    abs_b = sw
                #end
                #s = sw # uncomment this to start next lot with swept solution. Not reccomended since it usually gets you stuck in a local minimum.
                #println(cost(g,b,pun,avg))
                #println(cost(g,sw,pun,avg))
                #println("")
            #end
        end
       
        T = ϕ*T
        println(T)
    end

    v,abs_b
end

