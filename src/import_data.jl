export import_csv
"""
    import_csv(seq, file)
Import data from a CSV file. 

The 1st column should include the x data points (time or b-factor).
The 2nd column should be the y data points.
For imaginary numbers, the 2nd and 3rd columns can be used for 
the real and imaginary part.

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
without any unit conversions. 
Ensure that your x-axis is in seconds for relaxation, 
or in `seconds/metres² * 1e-9` for PFG experiments.
"""
function import_csv(t::Type{<:DataAxis}, file=pick_file(pwd()))

    x, y = import_csv(file)

    if isa(y, Vector{<:Complex})
        return autophase(ExperimentData(t(x),y))
    else
        return ExperimentData(t(x),y)
    end
end

function import_csv(file=pick_file(pwd()))
    data = readdlm(file, ',')

    if size(data,2) == 2

        x = vec(data[:, 1])
        y = vec(data[:, 2])
        return x, y

    elseif size(data,2) == 3

        x = vec(data[:, 1]) 
        y = complex.(vec(data[:, 2]) , vec(data[:, 3]))
        return x, y

    end
end

export import_dps
"""
    import_dps(seq, file)
Import data from a dps file, output from Bruker's minispec. 

The 1st column should be the data point index.
The 2nd column should include the x data point (time in ms).
The 3rd column should be the y data points.

The function reads the file and returns an `input1D` structure.
- `seq` is the 1D pulse sequence (e.g. IR, CPMG, PFG)
- `file` is the path to the dps file which contains the data (x, y) in two respective columns.

"""
function import_dps(t::Type{<:DataAxis}, file=pick_file(pwd()))

    data = readdlm(file, '\t')
    x = vec(data[:, 2]) .* 1e-3
    y = vec(data[:, 3])
    return ExperimentData(t(x), y)
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

export import_tnt, read_tnt_data, read_tnt_header
const TNT_MAGIC = r"^TNT1\.\d{3}$"

"""
Parse all TLV sections from a .tnt file.
Returns a Dict of tag => (offset, length, bool_flag, data).
"""
function parse_tnt_sections(filename::String)
    sections = Dict{String, NamedTuple}()
    open(filename, "r") do io
        # Validate magic number (8 bytes: e.g. "TNT1.000")
        magic = String(read(io, 8))
        occursin(TNT_MAGIC, strip(magic, '\0')) || error("Not a valid TNMR file: $magic")

        while !eof(io)
            # Need at least 12 bytes for a TLV header (tag=4, bool=4, len=4)
            pos_before = position(io)
            remaining = filesize(io) - pos_before
            remaining < 12 && break

            tag   = String(read(io, 4))
            bflag = read(io, UInt32)
            dlen  = read(io, UInt32)
            offset = position(io)

            if tag == "PSEQ"
                # PSEQ: the 4-byte dlen is NOT a length field — PSEQ has no length prefix.
                # The bytes we read as dlen are the start of the variable-length content.
                # Seek back 4 bytes and read everything to EOF.
                seek(io, offset - 4)
                data = read(io)
                sections[tag] = (offset = offset - 4, length = length(data),
                                 bool = Bool(bflag), data = data)
                break  # PSEQ is always last
            else
                data = read(io, dlen)
                sections[tag] = (offset = offset, length = Int(dlen),
                                 bool = Bool(bflag), data = data)
            end
        end
    end
    return sections
end

"""
    read_tnt_data(filename::String)

Read raw NMR data from a Tecmag .tnt file.
Returns a ComplexF32 array. Shape is determined by npts (all 4 dims),
which controls how much data was written, regardless of actual_npts.
"""
function read_tnt_data(filename::String)
    sections = parse_tnt_sections(filename)
    haskey(sections, "DATA") || error("No DATA section found")
    haskey(sections, "TMAG") || error("No TMAG section found")

    tmag = sections["TMAG"].data
    npts        = [reinterpret(Int32, tmag[4i-3 : 4i])[1]        for i in 1:4]
    actual_npts = [reinterpret(Int32, tmag[16 + 4i-3 : 16 + 4i])[1] for i in 1:4]

    data_bytes = sections["DATA"].data
    n_bytes    = length(data_bytes)
    n_complex  = n_bytes ÷ 8  # each ComplexF32 = 2 × Float32 = 8 bytes

    # Sanity check: data size must be exact
    n_bytes == n_complex * 8 || error("DATA size not a multiple of 8 bytes")

    # Build shape from actual_npts, treating 0/1 unused dims as 1
    # then drop trailing size-1 dims
    shape_full = Tuple(max(1, Int(n)) for n in actual_npts)
    prod(shape_full) == n_complex ||
        @warn "actual_npts product $(prod(shape_full)) ≠ data points $n_complex; reshaping to flat vector"

    data = copy(reinterpret(ComplexF32, data_bytes))

    if prod(shape_full) == n_complex
        # Drop trailing size-1 dimensions (e.g. (40000,1,1,1) → (40000,))
        active = findlast(n -> n > 1, collect(shape_full))
        shape  = isnothing(active) ? (n_complex,) : shape_full[1:active]
        return reshape(data, shape), actual_npts
    else
        # Fall back: return flat, let caller figure out shape
        return data, actual_npts
    end
end

"""
    read_tnt_header(filename::String)

Read TMAG parameters from a .tnt file and return them as a Dict.
Covers the full 1024-byte TECMAG structure per spec.
"""
function read_tnt_header(filename::String)
    sections = parse_tnt_sections(filename)
    haskey(sections, "TMAG") || error("No TMAG section found")
    tmag = sections["TMAG"].data
    length(tmag) >= 1024 || error("TMAG section too short: $(length(tmag)) bytes")

    io = IOBuffer(tmag)
    h  = Dict{String,Any}()

    # --- Points and scans (76 bytes) ---
    h["npts"]          = [read(io, Int32) for _ in 1:4]   # 16 bytes
    h["actual_npts"]   = [read(io, Int32) for _ in 1:4]   # 16 bytes
    h["acq_points"]    = read(io, Int32)                   #  4 bytes
    h["npts_start"]    = [read(io, Int32) for _ in 1:4]   # 16 bytes
    h["scans"]         = read(io, Int32)                   #  4 bytes
    h["actual_scans"]  = read(io, Int32)                   #  4 bytes
    h["dummy_scans"]   = read(io, Int32)                   #  4 bytes
    h["repeat_times"]  = read(io, Int32)                   #  4 bytes
    h["sadimension"]   = read(io, Int32)                   #  4 bytes
    h["samode"]        = read(io, Int32)                   #  4 bytes
    # space1[0] — zero bytes, nothing to skip
    # Running offset: 76

    # --- Field and frequencies (164 bytes cumulative, so 88 bytes here) ---
    h["magnet_field"]  = read(io, Float64)                 #  8
    h["ob_freq"]       = [read(io, Float64) for _ in 1:4]  # 32
    h["base_freq"]     = [read(io, Float64) for _ in 1:4]  # 32
    h["offset_freq"]   = [read(io, Float64) for _ in 1:4]  # 32
    h["ref_freq"]      = read(io, Float64)                 #  8
    h["NMR_frequency"] = read(io, Float64)                 #  8
    h["obs_channel"]   = read(io, Int16)                   #  2
    skip(io, 42)  # space2[42]
    # 8+32+32+32+8+8+2+42 = 164. ✓

    # --- Spectral width, dwell, filter (128 bytes) ---
    h["sw"]                 = [read(io, Float64) for _ in 1:4]  # 32
    h["dwell"]              = [read(io, Float64) for _ in 1:4]  # 32
    h["filter"]             = read(io, Float64)                  #  8
    h["experiment_time"]    = read(io, Float64)                  #  8
    h["acq_time"]           = read(io, Float64)                  #  8
    h["last_delay"]         = read(io, Float64)                  #  8
    h["spectrum_direction"] = read(io, Int16)                    #  2
    h["hardware_sideband"]  = read(io, Int16)                    #  2
    h["Taps"]               = read(io, Int16)                    #  2
    h["Type"]               = read(io, Int16)                    #  2
    h["bDigRec"]            = read(io, Int32)                    #  4  (BOOL)
    h["nDigitalCenter"]     = read(io, Int32)                    #  4
    skip(io, 16)  # space3[16]
    # 32+32+8+8+8+8+2+2+2+2+4+4+16 = 128 ✓

    # --- Hardware settings (20 bytes) ---
    h["transmitter_gain"]    = read(io, Int16)                   #  2
    h["receiver_gain"]       = read(io, Int16)                   #  2
    h["NumberOfReceivers"]   = read(io, Int16)                   #  2
    h["RG2"]                 = read(io, Int16)                   #  2
    h["receiver_phase"]      = read(io, Float64)                 #  8
    skip(io, 4)  # space4[4]
    # 2+2+2+2+8+4 = 20 ✓

    # --- Spinning (4 bytes) ---
    h["set_spin_rate"]       = read(io, UInt16)
    h["actual_spin_rate"]    = read(io, UInt16)

    # --- Lock (48 bytes) ---
    h["lock_field"]     = read(io, Int16)
    h["lock_power"]     = read(io, Int16)
    h["lock_gain"]      = read(io, Int16)
    h["lock_phase"]     = read(io, Int16)
    h["lock_freq_mhz"]  = read(io, Float64)
    h["lock_ppm"]       = read(io, Float64)
    h["H2O_freq_ref"]   = read(io, Float64)
    skip(io, 16)  # space5[16]
    # 2+2+2+2+8+8+8+16 = 48 ✓

    # --- VT (16 bytes) ---
    h["set_temperature"]    = read(io, Float64)
    h["actual_temperature"] = read(io, Float64)

    # --- Shim (88 bytes) ---
    h["shim_units"]  = read(io, Float64)
    h["shims"]       = [read(io, Int16) for _ in 1:36]
    h["shim_FWHM"]   = read(io, Float64)
    # 8 + 72 + 8 = 88 ✓

    # --- Bruker-specific (20 bytes) ---
    h["HH_dcpl_attn"]    = read(io, Int16)
    h["DF_DN"]           = read(io, Int16)
    h["F1_tran_mode"]    = [read(io, Int16) for _ in 1:7]
    h["dec_BW"]          = read(io, Int16)
    # 2+2+14+2 = 20 ✓

    # --- Gradient / misc (284 bytes) ---
    h["grd_orientation"] = String(read(io, 4))
    h["LatchLP"]         = read(io, Int32)
    h["grd_Theta"]       = read(io, Float64)
    h["grd_Phi"]         = read(io, Float64)
    skip(io, 264)  # space6[264]
    # 4+4+8+8+264 = 288  ← spec says 300 for this section including time vars below

    # --- Time variables (12 bytes, CTime/CTimeSpan = 4 bytes each) ---
    h["start_time"]   = read(io, UInt32)
    h["finish_time"]  = read(io, UInt32)
    h["elapsed_time"] = read(io, UInt32)
    # 4+4+8+8+264+4+4+4 = 300 ✓

    # --- Text variables (160 bytes) ---
    h["date"]         = rstrip(String(read(io, 32)), '\0')
    h["nucleus"]      = rstrip(String(read(io, 16)), '\0')
    h["nucleus_2D"]   = rstrip(String(read(io, 16)), '\0')
    h["nucleus_3D"]   = rstrip(String(read(io, 16)), '\0')
    h["nucleus_4D"]   = rstrip(String(read(io, 16)), '\0')
    h["sequence"]     = rstrip(String(read(io, 32)), '\0')
    h["lock_solvent"] = rstrip(String(read(io, 16)), '\0')
    h["lock_nucleus"] = rstrip(String(read(io, 16)), '\0')
    # 32+16+16+16+16+32+16+16 = 160 ✓
    # Total TECMAG = 76+164+128+20+4+48+16+88+20+300+160 = 1024 ✓

    return h
end

function read_tnt_echo_times(filename::String, len::Int, echotime::String)

    open(filename) do io
        readuntil(io, "DATA")
        read(io,Int32) # discard this (empty bytes)
        data_length = read(io,Int32) # how much data (in bytes)
        skip(io, data_length)
        while !eof(io)
            readuntil(io, echotime)
            if !isprint(peek(io, Char)) 
                tₑ = parse(Int, filter(isdigit , readuntil(io,"u"))) * 1e-6
                x = [ t * tₑ for t in 1:len ]
                return x
            end
        end
        throw("Could not find `$echotime` in the file.")
    end
end

function read_tnt_de0_times(filename, len::Int)

    open(filename) do io

        readuntil(io, "DATA")
        read(io,Int32) # discard this (empty bytes)
        data_length = read(io,Int32) # how much data (in bytes)
        skip(io, data_length)

        readuntil(io, "de0:2") # 1st instance is not relevant
        readuntil(io, "de0:2") # 2nd instance contains data
        read(io,Int32) # discard this (empty bytes)
        x = Vector{Float64}(undef,0)
        while !eof(io)
            try
                push!(x, (parse(Float64, readuntil(io, "m")) + 0.01) * 1e-3)
            catch
                break
            end
        end
        if length(x) > 1
            if length(x) >= len
                return x[1:len]
            else
                extra_points = len - length(x)
                extrapolated = x[end] .* ( (x[end] / x[end-1]) .^ (1:extra_points) )
                @warn "de0:2 points not enough, extrapolating logarithmically \
                to match size of the data."
                return [x ; extrapolated]
            end
        else
            throw("Could not read de0 times.")
        end
    end
end

"""
    import_tnt(seq, filename)
Read data from a tecmag .tnt file, and return a `input1D` or `input2D` object with 
all the relevant information.

This function is experimental, and since tecmag consoles work mainly with 
custom pulse sequences, it might not cover some cases. 
If it does not work for you, please submit an issue.

 Necessary (positional) arguments:
 - `seq` is the pulse sequence, i.e. `IR`, `SR`, `CPMG`, `(CPMG,IR)`, `(CPMG,SR)`.
- `filename` is the path of the file containing the data.

 Optional (keyword) arguments:
- `echotime` is the string that is used to denote the Echo time as
defined in the pulse sequence. Default value is "Echo_Time".

Calling this function without a filename argument, e.g.: `import_tnt(seq)`, 
will open a file dialogue to select the .tnt file.
"""
function import_tnt(
    seq::Union{Type, Tuple{Vararg{Type}}},
    filename::String =pick_file(pwd()); 
    echotime="Echo_Time"
    )

    header = read_tnt_header(filename)
    data, _ = read_tnt_data(filename)
    data = vec(sum(reshape(data , header["npts"][1], :), dims = 1))

    if seq isa Tuple 

        data = reshape(data , :, header["actual_npts"][2])

        options = ((CPMG,IR), (CPMG,SR))

        if seq in options

            x_dir = read_tnt_echo_times(filename, size(data,1), echotime)
            x_indir = read_tnt_de0_times(filename, size(data,2))

            return autophase(ExperimentData((seq[1](x_dir), seq[2](x_indir)), data))
        else
            error("Currently only $(options) are implemented here.")
        end

    else
        x = if seq in [IR,SR]
            read_tnt_de0_times(filename, length(data))
        elseif seq == CPMG
            read_tnt_echo_times(filename, length(data), echotime)

        else
            error("Currently only IR, SR and CPMG are impemented here.")
        end

        return autophase(ExperimentData(seq(x), data))
    end
end


export import_spinsolve
"""
    import_spinsolve(files)
Import data from a Spinsolve expert experiment. 
Two paths must be provided as follows (order is not important):
- `files` = [.../datafile.csv , .../acqu.par] 

Calling this function without an argument by typing `import_spinsolve()`
will open a file dialogue to select the files.

The function reads the acqu.par.bak file to get the 
acquisition parameters, and the .dat file to get the data. 
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

    if exp == "T1"
        return import_csv(IR, datafile)

    elseif exp == "T2"
        return import_csv(CPMG, datafile)

    elseif exp in ["PGSTE", "PGSE"]
        return import_csv(PFG, datafile)

    elseif exp == "T1IRT2"
        return spinsolve_read_IRCPMG(acqufile, datafile)

    elseif exp == "DT2"
        return spinsolve_read_PFGCPMG(acqufile, datafile)

    elseif exp == "T2T2"
        return spinsolve_read_CPMGCPMG(acqufile, datafile)

    else
        error("Unrecognised pulse sequence")
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

    return ExperimentData( (CPMG(t_direct), IR(t_indirect)), Data )
end

function spinsolve_read_PFGCPMG(acqufile, datafile)

    n_echoes = parse(Int32 ,read_acqu(acqufile, "nrEchoes"))
    t_echo = parse(Float64, read_acqu(acqufile, "echoTime")) * 1e-6

    ramp = parse(Float64, read_acqu(acqufile, "gradRamp")) * 1e-3
    Δ = parse(Float64, read_acqu(acqufile, "bDelta")) * 1e-3
    δ = ramp + parse(Float64, read_acqu(acqufile, "lDelta")) * 1e-3
    steps = parse(Int32, read_acqu(acqufile, "nrSteps"))
    gradmax = parse(Float64, read_acqu(acqufile, "gradMax")) * 1e-3

    g = range(0.001, gradmax ,steps) 

    γ = if read_acqu(acqufile, "nucleus") == "1H-1H"
        267.52218744e6; # (rad s^-1 T^-1) 
    else
        error("Unrecognised Nucleus in aqcu.par")
    end 

    t_direct = collect(1:n_echoes) * t_echo
    bfactor = g.^2 .* (γ^2 * δ^2 * (Δ - δ/3)) .* 1e-9

    Raw = readdlm(datafile, ' ')

    if size(Raw, 2) == 1
        Raw = readdlm(datafile, ',')
    end

    Data = collect(transpose(complex.(Raw[:, 1:2:end], Raw[:, 2:2:end])))

    return ExperimentData((CPMG(t_direct), PFG(bfactor)), Data)

end

function spinsolve_read_CPMGCPMG(acqufile, datafile)

    n_echoes = parse(Int32 ,read_acqu(acqufile, "nrEchoes"))
    t_echo = parse(Float64, read_acqu(acqufile, "echoTime")) * 1e-6
    τ_steps = parse(Int32, read_acqu(acqufile, "nrEchoSteps"))
    τ_min = parse(Float64, read_acqu(acqufile, "minEcho")) * 1e-3
    τ_max = parse(Float64, read_acqu(acqufile, "maxEcho")) * 1e-3

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

    # make it so that it's multiples of echotime
    #=map!(x -> round(x / t_echo) * t_echo, t_indirect)=#

    return ExperimentData( (CPMG(t_direct), CPMG(t_indirect)), Data )
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

    dimensions = [0, 0]
    open(filedir) do io

        readuntil(io, "TestType=")
        pulse_sequence_number = parse(Int16, readline(io))
        readuntil(io, "Dimensions=")
        dimensions .= parse.(Int16, split(readline(io), ','))
        readuntil(io, "[Data]")
        data = readdlm(io, '\t', Float64, skipstart=2)

        y_re = data[:, 3]
        y_im = data[:, 4]

        typedict = Dict(
            3 => "CPMG",
            7 => "IR",
            105 => "PFG",
            106 => "IRCPMG",
            108 => "PFGCPMG",
        )

        seq = typedict[pulse_sequence_number]

        if seq == "IRCPMG"
            return autophase(ExperimentData(
                (
                    CPMG(data[1:dimensions[1], 1] .* (1 / 1000) ),
                    IR(data[1:dimensions[1]:end, 2] .* (1 / 1000) ),
                ),
                reshape(complex.(y_re, y_im), dimensions[1], dimensions[2])
            ))
        elseif seq == "PFGCPMG"
            return autophase(ExperimentData(
                (
                    CPMG(data[1:dimensions[1], 1]),
                    PFG(data[1:dimensions[1]:end, 2] .* (1 / 1000)),
                ),
                reshape(complex.(y_re, y_im), dimensions[1], dimensions[2])
            ))
        elseif seq == "PFG"
            return autophase(ExperimentData(PFG(data[:, 1]), complex.(y_re, y_im)))
        elseif seq == "IR"
            return autophase(ExperimentData(IR(data[:, 1] .* (1 / 1000)), complex.(y_re, y_im)))
        elseif seq == "CPMG"
            return autophase(ExperimentData(CPMG(data[:, 1] .* (1 / 1000)), complex.(y_re, y_im)))
        end
    end
end
