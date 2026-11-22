# Measurement files: reading them while they are written, and recognising the committed
# prefix of an output file when a run is resumed from a checkpoint.

const FNV_OFFSET = 0xcbf29ce484222325
const FNV_PRIME = 0x00000100000001b3

"""
    fnv1a(bytes[, h]) -> UInt64

64-bit FNV-1a hash of `bytes` (a byte vector or a string), continued from `h`, so
hashing a file in pieces gives the same value as hashing it at once. Used to recognise
the committed prefix of an output file on resume: it detects accidental mismatches
(another run's file, a modified file), not deliberate ones.
"""
function fnv1a(bytes::AbstractVector{UInt8}, h::UInt64 = FNV_OFFSET)
    @inbounds for b in bytes
        h = (h ⊻ b) * FNV_PRIME
    end
    return h
end
fnv1a(s::AbstractString, h::UInt64 = FNV_OFFSET) = fnv1a(codeunits(s), h)

"""
    file_prefix_hash(path, n) -> UInt64

`fnv1a` of the first `n` bytes of the file at `path`; an error if the file is shorter.
"""
function file_prefix_hash(path::AbstractString, n::Integer)
    filesize(path) >= n || error("$path has $(filesize(path)) bytes, fewer than $n")
    h = FNV_OFFSET
    buf = Vector{UInt8}(undef, 1 << 20)
    open(path) do io
        left = n
        while left > 0
            k = readbytes!(io, buf, min(left, length(buf)))
            k > 0 || error("$path: unexpected end of file")
            h = fnv1a(view(buf, 1:k), h)
            left -= k
        end
    end
    return h
end

"""
    read_measurements(path) -> (header, params, X, dropped)

Read a measurement file of examples/glueball.jl or examples/blocked_glueball.jl: a
header line starting with `#`, then one row of numbers per configuration.

Only lines ending in a newline are read, so a file that is still being written can be
analysed: a final piece without a newline is a line the writer has not finished, and it
is ignored (`dropped` is then true). Every row that is read must parse completely as
numbers and have the width of the first row, and no second header may appear (two runs
joined together); otherwise an error names the line. Blank lines are skipped.

`header` is the header line, `params` its `key=value` pairs, and `X` the rows
(configurations × columns).
"""
function read_measurements(path::AbstractString)
    lines = split(read(path, String), '\n')
    dropped = !isempty(lines[end])          # no final newline: the last line is unfinished
    pop!(lines)
    (!isempty(lines) && startswith(lines[1], "#")) || error("$path: the first line is not a header")
    header = String(lines[1])
    params = Dict{String,String}()
    for p in split(header[2:end])
        occursin("=", p) || continue
        k, v = split(p, "="; limit = 2)
        params[String(k)] = String(v)
    end
    rows = Vector{Vector{Float64}}()
    width = 0
    for (n, line) in enumerate(lines)
        n == 1 && continue
        isempty(strip(line)) && continue
        startswith(line, "#") && error("$path line $n: a second header (two runs joined?)")
        tokens = split(line)
        row = Vector{Float64}(undef, length(tokens))
        for (j, t) in enumerate(tokens)
            x = tryparse(Float64, t)
            x === nothing && error("$path line $n, column $j: cannot parse \"$t\"")
            row[j] = x
        end
        width == 0 && (width = length(row))
        length(row) == width || error("$path line $n: $(length(row)) columns, expected $width")
        push!(rows, row)
    end
    isempty(rows) && error("$path: no complete data lines")
    return header, params, permutedims(reduce(hcat, rows)), dropped
end
