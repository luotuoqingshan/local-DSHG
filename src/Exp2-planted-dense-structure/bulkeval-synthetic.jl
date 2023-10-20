"""
    make_jobs(params, jobs)

Create jobs to be done in parallel.
Compute the anchored densest subhypergraph for each ϵ and R.
"""
function make_jobs(
    params::Vector{Tuple{Ti, Tf, Vector{Ti}, String}},
    jobs::RemoteChannel, 
) where {Ti <: Integer, Tf} 
    for t in params 
        put!(jobs, t)
    end
    for i = 1:length(workers())
        put!(jobs, (-1, -1, -1, ""))
    end
end


"""
    do_jobs(jobs, results, H, Ht)

Worker's job function, solving the anchored densest subhypergraph problem.
"""
function do_jobs(
    jobs::RemoteChannel,
    results::RemoteChannel,
    H::SparseMatrixCSC{Tf, Ti},
    Ht::SparseMatrixCSC{Tf, Ti},
) where {Ti <: Integer, Tf}
    while true
        # i-th cluster, j-th trial
        i, ϵ, R, penaltyType = take!(jobs)
        if i == -1
            break
        end
        type = "global"
        if ϵ >= 1.0
            type = "local"
        end
        if penaltyType in ["fracvol", "vol"]
            res = ADHSG_flow_density_improvement(H, Ht, R, ϵ, penaltyType;flowtol=1e-8, type=type)
            put!(results, (i, res))
        elseif penaltyType in ["WCE", "UCE"]
            res = ADSG_flow_density_improvement(H, Ht, R, ϵ;flowtol=1e-8, type=type, expansion=penaltyType)
            put!(results, (i, res))
        else
            error("unsupported penalty type.")
        end
    end
end


"""
    bulkeval_ADHSG(H, Ht, ϵ, cluster, Rs, penaltyType)

Solve the anchored densest subhypergraph problem for a bunch of Rs in parallel.
The cluster is the ground truth planted dense structure.
"""
function bulkeval_ADHSG(
    H::SparseMatrixCSC{Tf, Ti},
    Ht::SparseMatrixCSC{Tf, Ti},
    ϵ::Tf,
    cluster::Vector{Ti},
    Rs::Vector{Vector{Ti}},
    penaltyType::String,
) where {Ti <: Integer, Tf}
    ntrial = length(Rs)
    jobs = RemoteChannel(() -> Channel{Tuple}(ntrials + length(workers())))
    results = RemoteChannel(() -> Channel{Tuple}(ntrials))
    tasks = Tuple{Ti, Tf, Vector{Ti}, String}[]
    for j = 1:length(Rs)
        push!(tasks, (j, ϵ, Rs[j], penaltyType))
    end
    make_jobs(tasks, jobs)
    for p in workers()
        remote_do(do_jobs, p, jobs, results, H, Ht)
    end  
    njob = length(tasks)
    objs = zeros(Tf, ntrial)
    sizes = zeros(Ti, ntrial)
    Rprecisions = zeros(Tf, ntrial)
    F1scores = zeros(Tf, ntrial)
    dts = zeros(Tf, ntrial)
    while njob > 0
        i, res = take!(results)
        njob -= 1
        println("$njob left.")
        objs[i] = res["optval"]
        sizes[i] = length(res["optsol"])
        Rprecisions[i] = precision(res["optsol"], Rs[i]) 
        F1scores[i] = F1score(cluster, res["optsol"])
        dts[i] = res["total_dt"]
    end
    return Dict(
        "objs" => objs,
        "sizes" => sizes,
        "Rprecisions" => Rprecisions,
        "F1scores" => F1scores,
        "dts" => dts,
    )
end


"""
    plot_ADHSG!(axis, ϵs, result, label, col; yscale = identity)

Plot the F1 score for each ϵ.  Show the mean and standard error. 
"""
function plot_ADHSG!(
    axis,
    ϵs::Vector{Tf},
    result,
    label,
    col;
    yscale = identity,
) where{Tf}
    result_mean = mean(result, dims=2)[:, 1]
    lines!(axis, ϵs, result_mean, label=label, color=(col, 1), yscale=yscale)
    n = size(result, 1)
    lowcurve = zeros(Tf, n)
    upcurve = zeros(Tf, n)
    for i = 1:n
        stderr = sem(result[i, :])
        low = max(result_mean[i] - stderr, 0.0)
        up = min(result_mean[i] + stderr, 1.0)
        lowcurve[i] = low
        upcurve[i] = up
    end
    band!(axis, ϵs, upcurve, lowcurve, color = (col, 0.4), yscale=yscale)
    scatter!(axis, ϵs, result_mean, color=(col, 1), yscale=yscale)
end

