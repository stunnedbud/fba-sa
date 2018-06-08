include("metabolic_network.jl")

function calc_lot(m::Array{Float64,2}, B::Array{Float64,2}, sc::Array{Float64,1}, v::Array{Float64,2}, c::Array{Float64,2}, cons::Array{Float64,2}, T::Float64, L::Int, pun::Float64, pun2::Float64, pun3::Float64, plot_name::Int, step_size::Float64, precision::Float64)
    v1 = v
    best = v
    sc1 = sc
    cost_best = cost(v, m, c, cons, pun,pun2,pun3,precision)
    accepted_costs = [cost_best]
    count = 0
    attempts = 0
    max_attempts = L^2 #50L
    cost_v1 = cost(v1, m, c, cons, pun,pun2,pun3,precision)
    while count < L
        attempts += 1
        cost_v1 = cost(v1, m, c, cons, pun,pun2,pun3,precision)
        sc2, v2 = neighbor(sc1,v1,cons,B,step_size)
        cost_v2 = cost(v2, m, c, cons, pun,pun2,pun3,precision)
        if cost_v2 >= cost_v1 - T # accepts solution
            v1 = v2#linear_combination(B, sc2)#v2
            sc1 = sc2
            cost_v1 = cost_v2#cost(v1, m, c, cons, pun,pun2,pun3,precision)#cost_v2
            count += 1
            push!(accepted_costs, cost_v1)
        end
        if cost_v2 > cost_best
            best = v2
            cost_best = cost_v2
        end
        if attempts > max_attempts
            println("Exceeded max attempts")
            return sum(accepted_costs)/L, sc1, v1, best, true # true to end run
        end 
    end
    
    # Uncomment to save timeseries of accepted costs to file.
    # Deactivated it because my memory was running out.
    #out = ""
    #for c in accepted_costs
    #    out *= string(c)*"\n" 
    #end
    #open("results/plots/"*string(plot_name)*".txt","a") do f
    #    write(f, out)
    #end

    sum(accepted_costs)/L, sc1, v1, best, false
end


# "Main" function in this file. Calculates lot until average of accepted solutions improves, whereupon it
# decreases T by multiplying it by theta (cooling factor). Keeps doing that until T reaches epsilon.
# Returns last accepted solution and the absolute best found.
function acceptance_by_thresholds(m::Array{Float64,2}, B::Array{Float64,2}, sc::Array{Float64,1}, v::Array{Float64,2}, c::Array{Float64,2}, cons::Array{Float64,2}, T::Float64, L::Int, ε::Float64, ϕ::Float64, pun::Float64, pun2::Float64, pun3::Float64, plot_name::Int, step_size::Float64, precision::Float64)
    open("results/plots/"*string(plot_name)*".txt","w") do f
        write(f, "")
    end
    #srand(seed)
    p = 0
    abs_b = v # absolute best
    sc_best = sc
    b = v
    f = false # flag for exceeding max attempts
    while T > ε
        q = p - 1
        while p > q
            q = p
            p, sc, v, b, f = calc_lot(m, B, sc, v, c, cons, T, L, pun, pun2, pun3, plot_name, step_size, precision)    

            if cost(b,m,c,cons,pun,pun2,pun3,precision) > cost(abs_b,m,c,cons,pun,pun2,pun3,precision)
                abs_b = b
                sc_best = sc
            end
            
            if f # exceeded max attempts
                #println("run ending exceeded max attempts")
                return sc,v,sc_best,abs_b
            end
        end
       
        T = ϕ*T
        println(T)
    end

    sc,v,sc_best,abs_b
end

