export input1D
"""
    input1D(seq, x, y)
A structure containing the following elements:
- `seq` is the 1D pulse sequence (e.g. IR, CPMG, PGSE)
- `x`, the x values of the measurement (e.g time for relaxation or b-factor for diffusion).
- `y`, the y values of the measurement. 

"""
struct input1D
    seq::Type{<:pulse_sequence1D}
    x::AbstractVector{<:Real}
    y::AbstractVector
end


export input2D
"""
    input2D(seq, x, y)
A structure containing the following elements:
- `seq` is the 2D pulse sequence (e.g. IRCPMG)
- `x_direct` is the direct dimension acquisition parameter (e.g. the times when you aquire CPMG echoes).
- `x_indirect` is the indirect dimension acquisition parameter (e.g. all the delay times τ in your IR sequence).

"""
struct input2D
    seq::Type{<:pulse_sequence2D}
    x_direct::AbstractVector{<:Real}
    x_indirect::AbstractVector{<:Real}
    data::AbstractMatrix
end


export inv_out_1D
"""
    inv_out_1D(seq, x, y, xfit, yfit, X, f, r, SNR, α, wa)
A structure containing the following elements:
- `seq` is the 1D pulse sequence (e.g. IR, CPMG, PGSE)
- `x`, the x values of the measurement (e.g time for relaxation or b-factor for diffusion).
- `y`, the y values of the measurement.
- `xfit`, the x values of the fitted data.
- `yfit`, the y values of the fitted data.
- `X`, the x values of the inversion results.
- `f`, the inversion results.
- `r`, the residuals.
- `SNR`, the signal-to-noise ratio.
- `α`, the regularization parameter.
- `wa`, the weighted average of the inversion results (e.g. the mean relaxation time or diffusion coefficient).

"""
struct inv_out_1D
    seq::Type{<:pulse_sequence1D}
    x::Vector
    y::Vector
    xfit::Vector
    yfit::Vector
    X::Vector
    f::Vector
    r::Vector
    SNR::Real
    alpha::Real
    wa::Real
end



export inv_out_2D
"""
    inv_out_2D(seq, X_dir, X_indir, F, r, SNR, α, del_pol, sel_pol)
A structure containing the following elements:
- `seq` is the 2D pulse sequence (e.g. IRCPMG)
- `X_dir`, the x values of the direct dimension.
- `X_indir`, the x values of the indirect dimension.
- `F`, the inversion results as a matrix.
- `r`, the residuals.
- `SNR`, the signal-to-noise ratio.
- `alpha`, the regularization parameter.
- `filter`, apply a mask to filter out artefacts when plotting.
- `selections`, the selection masks (e.g. when you want to highlight some peaks in a T₁-T₂ map).

"""
struct inv_out_2D
    seq::Type{<:pulse_sequence2D}
    X_dir::Vector
    X_indir::Vector
    F::Matrix
    r::Vector
    SNR::Real
    alpha::Real
    filter::Matrix
    selections::Vector{Vector{Vector}}
end


export import_csv
"""
    import_csv(seq, file)
Import data from a CSV file.
The function reads the file and returns an `input1D` structure.
- `seq` is the 1D pulse sequence (e.g. IR, CPMG, PGSE)
- `file` is the path to the CSV file which contains the data (x, y) in two respective columns.
\
The function can be called without the seq argument, and the output will be the x and y vectors( `x,y =import_csv()`).
Alternatively, the function can also be called with only the seq argument, in which case a file dialog will open to select the file.
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
        p = readline(io)
    end
    return replace(p, "\"" => "")
end


export import_spinsolve
"""
    import_spinsolve(files)
Import data from a Spinsolve experiment. 
Two paths must be provided as follows (order is not important):
- `files` = [.../datafile.csv , .../acqu.par.bak] 
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
    elseif exp == "PGSTE"
        seq = PFG
    end

    if seq in [IR, CPMG, PFG]

        return import_csv(seq, datafile)

    elseif seq in [IRCPMG]

        # Read experiment parameters
        acqu = readdlm(acqufile)
        n_echoes = acqu[21, 3]
        t_echo = acqu[12, 3] * 1e-6
        τ_steps = acqu[36, 3]
        τ_min = acqu[20, 3] * 1e-3
        τ_max = acqu[19, 3] * 1e-3
        experiment = acqu[13, 3]

        Raw = readdlm(datafile, ' ')

        if size(Raw, 2) == 1
            Raw = readdlm(datafile, ',')
        end

        Data = collect(transpose(complex.(Raw[:, 1:2:end], Raw[:, 2:2:end])))

        ## Make time arrays
        # Time array in direct dimension
        t_direct = collect(1:n_echoes) * t_echo

        # Time array in direct dimension
        if acqu[18, 3] == "yes" # if log spacing is selected, do log array
            t_indirect = exp10.(range(log10(τ_min), log10(τ_max), τ_steps))
        else                   # otherwise, do a linear array
            t_indirect = collect(range(τ_min, τ_max, τ_steps))
        end

        return input2D(seq, t_direct, t_indirect, Data)

    end


end



function writeresults(results::Union{inv_out_1D,inv_out_2D}, dir)

    open(dir, "w") do io
        for field in fieldnames(typeof(results))
            data = getfield(results, field)
            datastring = isa(data, Array) ? join(data, ", ") : string(getfield(results, field))
            write(io, String(field) * " : " * datastring * "\n")
        end
    end

end

function readresults(dir::String)

    open(dir) do io

        readuntil(io, "seq : ")
        seq = eval(Meta.parse(readline(io)))
        supertype(seq) == pulse_sequence1D ? datatype = inv_out_1D : datatype = inv_out_2D

        for field in fieldnames(datatype)
            readuntil(io, field * " : ")

        end

    end

end



export import_geospec
"""
    import_geospec(dir)
Import data from a .txt format, as exported by Geospec instruments.

The function reads the relevant information, performs a phase correction on the data,
and returns an `input1D` or `input2D` structure.
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
