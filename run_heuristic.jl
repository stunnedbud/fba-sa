__precompile__()
module Run
include("simulated_annealing.jl")

# Main function. Reads all settings from given file, runs the simulated 
# annealing a specified number of times with random seeds.
# Outputs results to a console and to txt file (at results/[master-seed].txt)
function run(settings_file::String)

    # Read settings file, save into dictionary
    k = AbstractString[]    
    v = AbstractString[]    
    open(settings_file, "r") do f
        lines = readlines(f)
        for l in lines
            if startswith(l, "#") || (l=="\n") || length(l) < 2 # so we can comment out a line
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

    # Get basis for nullspace of S
    B = nullspace(S)
    
    # Get master seed
    master_seed = parse(Int, settings["seed"])
    srand(master_seed)
    
    # Get punishments
    pun = parse(Float64,settings["punishment"])
    pun2 = parse(Float64,settings["punishment2"])
    pun3 = parse(Float64,settings["punishment3"])
    
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
    sc_abs_best = Union{}
    best_run = -1
    precision =  parse(Float64,settings["balance_precision"])
    initial_range = parse(Float64,settings["initial_range"])
    initial_solution_range = parse(Float64,settings["initial_solution_range"])
    
    # Main loop
    for i in 1:num_runs
        current_seed =  10000seeds[i]^2 + 100seeds[i+num_runs]^2 + seeds[i+2num_runs]^2 + master_seed 
        if num_runs == 1 # this way we detect when we want to replicate a given seed's results
            #println()
            current_seed = master_seed # we use the given seed directly when settings ask for a single run
        end
        srand(current_seed)
        #current_seed = master_seed
        region = rand(0:initial_solution_range)
        if rand() > 0.5 
            region = -region
        end
        sc,v = stable_solution(B, region, initial_range)
        if abs_best == Union{}
            abs_best = v # to initialize it in something not null
        end
        
        out = "\n\nRun #"*string(i)*"\n"
        out *= "Seed: "
        out *= string(current_seed)*"\n"
        out *= "Initial flux vector ["*string(feasible(v, S, cons, precision))*", "*string(count_constraints_fails(v,cons,precision))*" constraints failed] ("*string(cost(v, S, c, cons, pun, pun2, pun3, precision))*"): "
        out *= "\nb: "
        for j in sc
            out *= string(j)*","
        end
        out = out[1:end-1]*"\nv: "
        for j in v
            out *= string(j)*","
        end
        out = out[1:end-1]*"\n\n"
    
        # Run heuristic
        sc_last, last, sc_best, best = acceptance_by_thresholds(S, B, sc, v, c, cons, parse(Float64,settings["T"]), parse(Int,settings["L"]), parse(Float64,settings["epsilon"]), parse(Float64,settings["theta"]), pun, pun2, pun3, current_seed, parse(Float64,settings["step_size"]), precision)
        
        # Update absolute best found
        if cost(best, S, c, cons, pun, pun2, pun3, precision) > cost(abs_best, S, c, cons, pun, pun2, pun3, precision)
            abs_best = best
            sc_abs_best = sc_best
            best_run = i
        end
        
        out *= "Last ["*string(feasible(last, S, cons, precision))*", "*string(count_constraints_fails(last,cons,precision))*" constraints failed] ("*string(cost(last, S, c, cons, pun, pun2, pun3, precision))*"):"
        out *= "\nb: "
        for j in sc_last
            out *= string(j)*","
        end

        out = out[1:end-1]*"\nv: "
        for j in last
            out *= string(j)*","
        end
        out = out[1:end-1]*"\n" # removes last comma
        
        out *= "Best ["*string(feasible(best, S, cons, precision))*", "*string(count_constraints_fails(best,cons,precision))*" constraints failed] ("*string(cost(best, S, c, cons, pun,pun2, pun3, precision))*"):"
        out *= "\nb: "
        for j in sc_best
            out *= string(j)*","
        end

        out = out[1:end-1]*"\nv: "
        for j in best
            out *= string(j)*","
        end
        out = out[1:end-1]*"\n" # removes last comma

        print(out)
        
        # Append run results to file
        open("results/"*string(master_seed)*".txt", "a") do f    
            write(f, out)
        end
    end
    
    # Final results
    out = "\n############################################################\n"
    out *= "Best solution found: \nb:"
    for j in sc_abs_best
        out *= string(j)*","
    end
    out = out[1:end-1]*"\nv:" # removes last comma
    for j in abs_best
        out *= string(j)*","
    end
    out = out[1:end-1]*"\n"
    out *= "With biomass "*string(dot(abs_best[:], c[:]))*"\n"
    out *= "Feasible? ["*string(feasible(abs_best, S, cons, precision))*", "*string(count_constraints_fails(abs_best,cons,precision))*" constraints failed]\n"
    out *= "From run #"*string(best_run)*"\n"

    open("results/"*string(master_seed)*".txt","a") do f
        write(f, out)
    end
    println(out)
end
end#module
