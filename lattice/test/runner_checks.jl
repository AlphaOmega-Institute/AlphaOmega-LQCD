# Checks of examples/blocked_glueball.jl as a program: interrupted and resumed chains,
# the refusals that protect checkpointed measurements, and reading the output while it is
# being written.
#
#   julia --project=lattice lattice/test/runner_checks.jl
#
# Small chains (8³ × 8 fine lattices, two stout levels, two α) run in subprocesses; the
# whole script takes a few minutes, most of it Julia start-up.

using Test, Serialization, LatticeYM

const ROOT = dirname(dirname(@__DIR__))
const RUNNER = joinpath(ROOT, "lattice", "examples", "blocked_glueball.jl")
const SMALL = ["--L", "8", "--T", "8", "--therm", "5", "--levels", "1,2", "--alphas", "0.0,0.3"]
const DIR = mktempdir()

file(name) = joinpath(DIR, name)
cmd(args...) = `$(Base.julia_cmd()) --project=$(joinpath(ROOT, "lattice")) -t 2 $RUNNER $SMALL $(collect(args))`
# run to completion; return (succeeded, standard error)
function runner(args...)
    err = file("stderr.txt")
    ok = success(pipeline(cmd(args...); stdout = devnull, stderr = err))
    return ok, read(err, String)
end
done_in(ckpt) = isfile(ckpt) ? deserialize(ckpt).done : 0
nlines(path) = isfile(path) ? count(==(UInt8('\n')), read(path)) : 0
# start a run, wait until `ready()` holds and then `delay` seconds more, and kill it hard. The
# output is flushed only at checkpoints, so the configurations since the last one are lost, as
# in a crash.
function run_and_kill(ready, delay, args...)
    p = run(pipeline(cmd(args...); stdout = devnull, stderr = devnull); wait = false)
    t = time()
    while !ready()
        process_exited(p) && error("the run finished before it could be interrupted")
        time() - t < 600 || (kill(p); error("timed out"))
        sleep(0.05)
    end
    sleep(delay)
    process_exited(p) && error("the run finished before it could be interrupted")
    kill(p, Base.SIGKILL)
    wait(p)
end

@testset "blocked_glueball.jl: checkpoints, resume, refusals, live reading" begin
    ref, run1, ck1 = file("ref.txt"), file("run.txt"), file("run.ckpt")

    # an uninterrupted chain of 75 configurations is the reference; while it runs, every read of
    # the growing file returns only complete rows, equal to rows of the finished file
    p = run(pipeline(cmd("--confs", "75", "--out", ref); stdout = devnull, stderr = devnull); wait = false)
    snapshots = Matrix{Float64}[]
    while !process_exited(p)
        if nlines(ref) >= 2
            _, _, X, _ = read_measurements(ref)
            push!(snapshots, X)
        end
        sleep(0.02)
    end
    @test success(p)
    _, _, Xref, dropped = read_measurements(ref)
    @test size(Xref, 1) == 75 && !dropped
    @test length(snapshots) > 10
    @test all(X -> X == Xref[1:size(X, 1), :], snapshots)

    # killed after the checkpoint at 25, resumed, killed after the one at 50, resumed again:
    # the output is byte-identical to the uninterrupted chain
    run_and_kill(() -> done_in(ck1) == 25, 1.0, "--confs", "75", "--out", run1, "--checkpoint", ck1)
    ok, err = runner("--resume", ck1, "--out", run1, "--confs", "40")         # a shorter target first
    @test ok && occursin("resuming after configuration 25", err)
    @test done_in(ck1) == 40                                                 # final checkpoint, not a multiple of 25
    run_and_kill(() -> done_in(ck1) == 50, 1.0, "--resume", ck1, "--out", run1, "--confs", "75")
    ok, err = runner("--resume", ck1, "--out", run1, "--confs", "75")
    @test ok && occursin("resuming after configuration 50", err)
    @test read(run1) == read(ref)

    # a fresh run does not overwrite an existing output or checkpoint file
    before = read(run1)
    ok, err = runner("--confs", "75", "--out", run1, "--checkpoint", file("new.ckpt"))
    @test !ok && occursin("refusing to overwrite", err) && read(run1) == before
    ok, err = runner("--confs", "60", "--out", file("new.txt"), "--checkpoint", ck1)
    @test !ok && occursin("refusing to overwrite", err) && !isfile(file("new.txt"))
    ok, err = runner("--confs", "60", "--checkpoint", file("new.ckpt"))
    @test !ok && occursin("need --out", err)

    # a checkpoint of one chain cannot be resumed onto the output of another (seed 10 vs seed 11)
    o10, c10, o11, c11 = file("s10.txt"), file("s10.ckpt"), file("s11.txt"), file("s11.ckpt")
    @test runner("--seed", "10", "--confs", "25", "--out", o10, "--checkpoint", c10)[1]
    @test runner("--seed", "11", "--confs", "25", "--out", o11, "--checkpoint", c11)[1]
    b11 = read(o11)
    ok, err = runner("--resume", c10, "--out", o11, "--confs", "50")
    @test !ok && occursin("different header", err) && read(o11) == b11
    # nor onto a file with the right header whose committed rows were changed, or cut short
    lines = readlines(o10)
    lines[3] = replace(lines[3], r"\d(?=e)" => s -> string((parse(Int, s) + 1) % 10); count = 1)
    write(o11, join(lines, "\n") * "\n")
    ok, err = runner("--resume", c10, "--out", o11, "--confs", "50")
    @test !ok && occursin("differ from those recorded", err)
    write(o11, read(o10)[1:end-100])
    ok, err = runner("--resume", c10, "--out", o11, "--confs", "50")
    @test !ok && occursin("fewer than", err)
    # physics options given on resume must agree with the checkpoint
    ok, err = runner("--resume", c10, "--out", o10, "--confs", "50", "--seed", "11")
    @test !ok && occursin("differs from the checkpoint", err)

    # lines written after the checkpoint, including an unfinished one, are cut off in place;
    # the resumed chain (physics taken from the checkpoint) matches an uninterrupted one
    open(o10, "a") do f
        write(f, "0.59 1.234e-01 5.6e-02\n0.58 9.9e-")
    end
    ok, err = runner("--resume", c10, "--out", o10, "--confs", "40")
    @test ok
    ok2, _ = runner("--seed", "10", "--confs", "40", "--out", file("s10ref.txt"))
    @test ok2 && read(o10) == read(file("s10ref.txt"))
end
