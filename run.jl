__precompile__()
module Run
include("simulated_annealing.jl")

# Main function. Reads all settings from given file, runs the simulated 
# annealing a specified number of times with random seeds.
# Outputs results to a console and to txt file (at results/[master-seed].txt)
function run(settings_file::UTF8String)

    # Read settings file, save into dictionary
    k = AbstractString[]    
    v = AbstractString[]    
    open(settings_file, "r") do f
        lines = readlines(f)
        for l in lines
            if startswith(l, "#") | (l=="\n") # so we can comment out a line
                continue
            end
            str = split(strip(l), "=")
            push!(k, str[1])
            push!(v, str[2])
        end
    end
    settings = Dict(zip(k,v))

    println("Settings")
    for (k,v) in settings
        println(k*": "*string(v))
    end
    
    # Create needed dirs for results if they don't exist
    if !isdir("results")
        mkdir("results")
    end
    if !isdir("results/plots")
        mkdir("results/plots")
    end
    
    # Parse matrix from csv
    S, reactions, metabolites  = parse_matrix(settings["network_source_file"])
    # Parse objective function
    c = parse_objective(settings["objective_source_file"])
    # Parse constraints
    cons = parse_constraints(settings["constraints_source_file"])
    
    # Get master seed
    master_seed = parse(Int, settings["seed"])
    srand(master_seed)
    
    # Get punishment
    pun = parse(Float64,settings["punishment"])
    
    # Results file
    out = "RUN RESULTS\n\n#######################################\nInitial settings\n"
    for (k,v) in settings
        out *= k*": "*string(v)*"\n" 
    end 
    out *= "#######################################\n\n"
    open("results/"*string(master_seed)*".txt","w") do f
        write(f, out)
    end
    
    # Number of runs
    num_runs = parse(Int, settings["runs"])

    # Generate seeds
    seeds = randperm(num_runs*3)

    # Run heuristic
    abs_best = Union{}
    best_run = -1
    
    # Main loop
    for i in 1:num_runs
        current_seed =  10000seeds[i]^2 + 100seeds[i+num_runs]^2 + seeds[i+2num_runs]^2 + master_seed
        v = random_solution(cons)
        if abs_best == Union{}
            abs_best = v # to initialize it in something not null
        end
        
        out = "Run #"*string(i)*"\n"
        out *= "Seed: "
        out *= string(current_seed)*"\n"
        out *= "Initial flux vector ["*string(feasible(v, S, cons))*", "*string(stability_sum(v,S))*"] ("*string(cost(v, S, c, cons, pun))*"): "
        for j in v
            out *= string(j)*","
        end
        out = out[1:end-1]*"\n\n"
    
        # Run heuristic
        last, best = acceptance_by_thresholds(S, v, c, cons, parse(Float64,utf8(settings["T"])), parse(Int,utf8(settings["L"])), parse(Float64,utf8(settings["epsilon"])), parse(Float64,utf8(settings["theta"])), pun, current_seed)
        
        # Update absolute best found
        if cost(best, S, c, cons, pun) > cost(abs_best, S, c, cons, pun)
            abs_best = best
            best_run = i
        end
        
        out *= "Last ["*string(feasible(last, S, cons))*", "*string(stability_sum(last,S))*"] ("*string(cost(last, S, c, cons, pun))*"): "
        for j in last
            out *= string(j)*","
        end
        out = out[1:end-1]*"\n" # removes last comma
        out *= "Best ["*string(feasible(best, S, cons))*", "*string(stability_sum(best,S))*"] ("*string(cost(best, S, c, cons, pun))*"): "
        for j in best
            out *= string(j)*","
        end
        out = out[1:end-1]*"\n\n" # remove last comma
        print(out)
        
        # Append run results to file
        open("results/"*string(master_seed)*".txt", "a") do f    
            write(f, out)
        end
    end
    
    # Final results
    out = "\n############################################################\n"
    out *= "Best solution found: \n"
    for j in abs_best
        out *= string(j)*","
    end
    out = out[1:end-1]*"\n" # removes last comma
    out *= "With biomass "*string(cost(abs_best, S, c, cons, pun))*"\n"
    out *= "Feasible? ["*string(feasible(abs_best, S, cons))*", "*string(stability_sum(abs_best,S))*"]\n"
    out *= "From run #"*string(best_run)*"\n"

    open("results/"*string(master_seed)*".txt","a") do f
        write(f, out)
    end
    println(out)
end
end#module
