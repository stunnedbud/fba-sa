
# Reads a csv file containing stoichiometric matrix for a metabolic network
# Returns parsed matrix, as well as one list of reaction ids and one for metabolite ids.
function parse_matrix(csv_source_file::AbstractString)   
    network = None
    reactions = AbstractString[]
    metabolites = AbstractString[]
     
    open(csv_source_file) do f
        lines = readlines(f)
        firstline = true
        i = 1
        for l in lines
            if firstline # first line contains column headers (ie. reaction ids)
                words = split(l, ",")
                words = map(utf8, words)
                for w in words
                    if firstline # skips first word
                        firstline = false
                        continue
                    end
                    push!(reactions, w) 
                end
                network = fill(0.0, (length(lines)-1, length(reactions))) 
                continue
            end
            
            mname = true
            words = split(l, ",")
            j = 1
            for w in words
                if mname # first column contains metabolites
                    push!(metabolites,w)
                    mname = false
                    continue
                end
                network[i, j] = parse(Float64, w)
                j+=1
            end
            i+=1
        end
    end
    
    network, reactions, metabolites
end

# Parses objective weights from file. Returns the c desired for the system:
# maximize v (dot) c    s.t.     S * v = 0
# Where S is the stoichiometry matrix and v is the flux vector we want to maximize.
# c represents the objective function, which transforms v into a biomass number.
function parse_objective(csv_source_file::AbstractString)
    c = None    
    open (csv_source_file) do f
        lines = readlines(f)
        for l in lines
            words = split(l, ",")
            c = fill(0.0, (length(words), 1))
            i = 1
            for w in words
                c[i] = parse(Float64, w)
                i+=1
            end
        end
    end
    return c
end

# Returns a matrix with minimum and maximum value allowed for each metabolite.
# Each row corresponds to values for each metabolite; the first column has the
# minimums and the second the maximums.
function parse_constraints(csv_source_file::AbstractString)
    c = None    
    open (csv_source_file) do f
        lines = readlines(f)
        c = fill(0.0, (length(lines)-1, 2))
        firstline = true
        i = 1
        for l in lines
            if firstline
                firstline = false
                continue
            end
            
            words = split(l, ",")
            c[i][1] = parse(Float64, words[2])
            c[i][2] = parse(Float64, words[3])
            i+=1
        end
    end
    return c
end


