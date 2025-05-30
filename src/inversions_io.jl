export import_csv
"""
    import_csv(seq, file)
Import data from a CSV file.
The function reads the file and returns an `input1D` structure.
- `seq` is the 1D pulse sequence (e.g. IR, CPMG, PFG)
- `file` is the path to the CSV file which contains the data (x, y) in two respective columns.


The function can be called without the seq argument,
and the output will be the x and y vectors 
(use it like, `x,y =import_csv()`).
Alternatively, the function can also be called with only the seq argument,
in which case a file dialog will open to select the file
(use it like, `data = import_csv(IR)`)


Please note that this function will just import the data as is,
without any unit conversions. Ensure that your x-axis is in SI.
"""
function import_csv(seq::Type{<:pulse_sequence1D}, file=pick_file(pwd()))
    x, y = import_csv(file)
    return input1D(seq, x, y)
end

function import_csv(file=pick_file(pwd()))
    data = readdlm(file, ',')
    x = vec(data[:, 1])
    y = vec(data[:, 2])
    return x, y
end


function read_acqu(filename, parameter)

    p = ""
    open(filename) do io

        readuntil(io, parameter * " = ")

        if !eof(io)
            p = readline(io)
        else
            seekstart(io)
            readuntil(io, parameter * "= ")
            p = readline(io)
            
        end
    end

    if p == ""
        throw("There seems to be an issue in reading the acqu file")
    end

    return replace(p, "\"" => "")
end


export import_spinsolve
"""
    import_spinsolve(files)
Import data from a Spinsolve experiment. 
Two paths must be provided as follows (order is not important):
- `files` = [.../datafile.csv , .../acqu.par] 
\
Calling this function without an argument by typing `import_spinsolve()` will open a file dialog to select the files.
The function reads the acqu.par.bak file to get the acquisition parameters, and the .dat file to get the data. 
The function returns an `input2D` structure.
"""
function import_spinsolve(files=pick_multi_file(pwd()))

    if length(files) == 1
        error("Please select both the .dat and 'acqu.par' files.")
    elseif length(files) > 2
        error("Please select only two files (data file and 'acqu.par'")
    end


    acqufile = files[contains.(files, "acqu")][1]
    datafile = files[.!contains.(files, "acqu")][1]

    exp = read_acqu(acqufile, "experiment")

    if exp == "T1IRT2"
        seq = IRCPMG
    elseif exp == "T1"
        seq = IR
    elseif exp == "T2"
        seq = CPMG
    elseif exp in ["PGSTE", "PGSE"]
        seq = PFG
    elseif exp == "DT2"
        seq = PFGCPMG
    else
        error("Unrecognised pulse sequence")
    end

    if seq in [IR, CPMG, PFG]

        return import_csv(seq, datafile)

    elseif seq in [IRCPMG]

        return spinsolve_read_IRCPMG(acqufile, datafile)

    elseif seq in [PFGCPMG]

        return spinsolve_read_PFGCPMG(acqufile, datafile)
    end

end

function spinsolve_read_IRCPMG(acqufile, datafile)
    
    n_echoes = parse(Int32 ,read_acqu(acqufile, "nrEchoes"))
    t_echo = parse(Float64, read_acqu(acqufile, "echoTime")) * 1e-6
    τ_steps = parse(Int32, read_acqu(acqufile, "tauSteps"))
    τ_min = parse(Float64, read_acqu(acqufile, "minTau")) * 1e-3
    τ_max = parse(Float64, read_acqu(acqufile, "maxTau")) * 1e-3

    Raw = readdlm(datafile, ' ')

    if size(Raw, 2) == 1
        Raw = readdlm(datafile, ',')
    end

    Data = collect(transpose(complex.(Raw[:, 1:2:end], Raw[:, 2:2:end])))

    ## Make time arrays
    # Time array in direct dimension
    t_direct = collect(1:n_echoes) * t_echo

    # Time array in direct dimension
    if  read_acqu(acqufile, "logspace") == "yes" # if log spacing is selected, do log array

        t_indirect = exp10.(range(log10(τ_min), log10(τ_max), τ_steps))
    else                   # otherwise, do a linear array

        t_indirect = collect(range(τ_min, τ_max, τ_steps))
    end

    return input2D(IRCPMG, t_direct, t_indirect, Data)
end

function spinsolve_read_PFGCPMG(acqufile, datafile)
    
    n_echoes = parse(Int32 ,read_acqu(acqufile, "nrEchoes"))
    t_echo = parse(Float64, read_acqu(acqufile, "echoTime")) * 1e-6
    Δ = parse(Int32, read_acqu(acqufile, "bDelta")) * 1e-3
    δ = parse(Int32, read_acqu(acqufile, "lDelta")) * 1e-3
    steps = parse(Int32, read_acqu(acqufile, "nrSteps"))
    gradmax = parse(Int32, read_acqu(acqufile, "gradMax")) * 1e-3
    g = range(0,steps,steps+1)[2:end] * (gradmax / steps) 

    γ = if read_acqu(acqufile, "nucleus") == "1H-1H"
        267.52218744e6; # (rad s^-1 T^-1) 
    else
        error("Unrecognised Nucleus in aqcu.par")
    end 

    t_direct = collect(1:n_echoes) * t_echo
    bfactor = g.^2 .* (γ^2 * δ^2 * (Δ - δ/3))

    Raw = readdlm(datafile, ' ')

    if size(Raw, 2) == 1
        Raw = readdlm(datafile, ',')
    end

    Data = collect(transpose(complex.(Raw[:, 1:2:end], Raw[:, 2:2:end])))

    return input2D(PFGCPMG, t_direct, bfactor, Data)

end


export import_geospec
"""
    import_geospec(dir)
Import data from a .txt format, as exported by Geospec instruments.

The function reads the relevant information, performs a phase correction on the data,
and returns an `input1D` or `input2D` structure.

Calling this function without an argument by typing `import_geospec()` 
will open a file dialog to select the .txt file.
"""
function import_geospec(filedir::String=pick_file(pwd()))

    # cd(dirname(filedir))

    data = []
    pulse_sequence_number::Int16 = 0
    dimensions = [0, 0]

    open(filedir) do io

        readuntil(io, "TestType=")
        pulse_sequence_number = parse(Int16, readline(io))
        readuntil(io, "Dimensions=")
        dimensions .= parse.(Int16, split(readline(io), ','))
        readuntil(io, "[Data]")
        data = readdlm(io, '\t', Float64, skipstart=2)
    end

    y_re = data[:, 3]
    y_im = data[:, 4]

    typedict = Dict(
        3 => CPMG,
        7 => IR,
        105 => PFG,
        106 => IRCPMG,
        108 => PFGCPMG
    )

    seq = typedict[pulse_sequence_number]

    if seq in [IR, IRCPMG]
        y_re, y_im, ϕ = autophase(y_re, y_im, -1)
    else
        y_re, y_im, ϕ = autophase(y_re, y_im, 1)
    end

    display("Data phase corrected by $(round(ϕ,digits=3)) radians.")

    if seq == IRCPMG

        return input2D(IRCPMG, data[1:dimensions[1], 1] .* (1 / 1000), data[1:dimensions[1]:end, 2] .* (1 / 1000), reshape(complex.(y_re, y_im), dimensions[1], dimensions[2]))

    elseif seq == PFGCPMG

        return input2D(PFGCPMG, data[1:dimensions[1], 1], data[1:dimensions[1]:end, 2] .* (1 / 1000), reshape(complex.(y_re, y_im), dimensions[1], dimensions[2]))

    elseif seq == PFG

        return input1D(seq, data[:, 1], complex.(y_re, y_im))

    elseif seq in [IR, CPMG]

        return input1D(seq, data[:, 1] .* (1 / 1000), complex.(y_re, y_im)) # Converts time to seconds
    end
end
