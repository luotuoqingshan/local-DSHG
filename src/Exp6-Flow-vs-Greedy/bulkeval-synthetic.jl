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
            flow_res = ADHSG_flow_density_improvement(H, Ht, R, ϵ, penaltyType;flowtol=1e-8, type=type)
            _, n = size(H)
            p = frac_volume_penalty(H, Ht, R, n, ϵ)
            greedy_res = greedy_peeling(H, Ht, p) 
            put!(results, (i, flow_res, greedy_res))
        else
            error("For this experiment, we use fracvol/vol penalty only.")
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
    flow_objs = zeros(Tf, ntrial)
    flow_sizes = zeros(Ti, ntrial)
    flow_Rprecisions = zeros(Tf, ntrial)
    flow_F1scores = zeros(Tf, ntrial)
    flow_dts = zeros(Tf, ntrial)
    greedy_objs = zeros(Tf, ntrial)
    greedy_sizes = zeros(Ti, ntrial)
    greedy_Rprecisions = zeros(Tf, ntrial)
    greedy_F1scores = zeros(Tf, ntrial)
    greedy_dts = zeros(Tf, ntrial)
    while njob > 0
        i, flow_res, greedy_res = take!(results)
        njob -= 1
        println("$njob left.")
        flow_objs[i] = flow_res["optval"]
        flow_sizes[i] = length(flow_res["optsol"])
        flow_Rprecisions[i] = precision(flow_res["optsol"], Rs[i]) 
        flow_F1scores[i] = F1score(cluster, flow_res["optsol"])
        flow_dts[i] = flow_res["total_dt"]
        greedy_objs[i] = greedy_res["optval"]
        greedy_sizes[i] = length(greedy_res["optsol"])
        greedy_Rprecisions[i] = precision(greedy_res["optsol"], Rs[i]) 
        greedy_F1scores[i] = F1score(cluster, greedy_res["optsol"])
        greedy_dts[i] = greedy_res["total_dt"]
    end
    return Dict(
        "flow_objs" => flow_objs,
        "flow_sizes" => flow_sizes,
        "flow_Rprecisions" => flow_Rprecisions,
        "flow_F1scores" => flow_F1scores,
        "flow_dts" => flow_dts,
        "greedy_objs" => greedy_objs,
        "greedy_sizes" => greedy_sizes,
        "greedy_Rprecisions" => greedy_Rprecisions,
        "greedy_F1scores" => greedy_F1scores,
        "greedy_dts" => greedy_dts,
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

