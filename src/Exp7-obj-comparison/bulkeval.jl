"""
    make_jobs(params, jobs)

Put (i, j, ϵ, R) tuples into the jobs RemoteChannel.
"""
function make_jobs(
    params::Vector{Tuple{Ti, Ti, Tf, Vector{Ti}, String}},
    jobs::RemoteChannel, 
) where {Ti <: Integer, Tf} 
    for t in params 
        put!(jobs, t)
    end
    for i = 1:length(workers())
        put!(jobs, (-1, -1, -1, -1, ""))
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
        # the i-th ϵ, j-th R
        i, j, ϵ, R, penaltyType = take!(jobs)
        # terminate when all jobs are done
        if i == -1
            break
        end
        type = "global"
        if ϵ >= 1.0
            type = "local"
        end
        if penaltyType in ["fracvol", "vol"]
            res = ADHSG_flow_density_improvement(H, Ht, R, ϵ, penaltyType;flowtol=1e-8, type=type)
            put!(results, (i, j, res)) 
        elseif penaltyType in ["WCE", "UCE"]
            res = ADSG_flow_density_improvement(H, Ht, R, ϵ;flowtol=1e-8, type=type, expansion=penaltyType)
            put!(results, (i, j, res)) 
        else
            error("unsupported penalty type.")
        end
    end
end


"""
    bulkeval_ADHSG(H, Ht, Rs, ϵs, penaltyType; [dataset="HSBM", savepath=homedir()*"/local-DHSG/results", filename="placeholder"]

Function for solving the Anchored Densest HyperSubgraph Problem 
for a bunch of ϵs and Rs in parallel.
"""
function bulkeval_ADHSG(
    H::SparseMatrixCSC{Tf, Ti},
    Ht::SparseMatrixCSC{Tf, Ti},
    Rs::Any,
    ϵs::Vector{Tf},
    penaltyType::String;
    dataset::String="HSBM",
    savepath::String=homedir()*"/local-DHSG/results",
    filename::String="placeholder",
) where {Ti <: Integer, Tf}
    ntrials = length(Rs)
    nϵ = length(ϵs)

    # set up jobs and results remote channel
    jobs = RemoteChannel(() -> Channel{Tuple}(ntrials * nϵ + length(workers())))
    results = RemoteChannel(() -> Channel{Tuple}(ntrials * nϵ))
    tasks = Tuple{Ti, Ti, Tf, Vector{Ti}, String}[]
    for i = 1:nϵ
        for j = 1:ntrials
            push!(tasks, (i, j, ϵs[i], Rs[j], penaltyType))
        end
    end
    make_jobs(tasks, jobs)
    for p in workers()
        remote_do(do_jobs, p, jobs, results, H, Ht)
    end  
    njob = length(tasks)
    Volpenalizeddensities = zeros(Tf, nϵ, ntrials)
    FracVolpenalizeddensities = zeros(Tf, nϵ, ntrials)
    while njob > 0
        i, j, res = take!(results)
        njob -= 1
        println("$njob left.")
        order = H_order(Ht)
        vol_penalty = volume_penalty(H, Rs[j], size(H, 2), ϵs[i])
        fracvol_penalty = frac_volume_penalty(H, Ht, Rs[j], size(H, 2), ϵs[i])
        Volpenalizeddensities[i, j] = edensity(H, order, vol_penalty, res["optsol"])
        FracVolpenalizeddensities[i, j] = edensity(H, order, fracvol_penalty, res["optsol"])
    end
    res = Dict(
        "epss" => ϵs,
        "Rs" => Rs,
        "ntrial" => length(Rs),
        "Volpenalizeddensities" => Volpenalizeddensities,
        "FractionalVolpenalizeddensities" => FracVolpenalizeddensities,
    )
    # save the results
    matwrite(savepath*"/"*dataset*"/"*filename*".mat", res)
    return res
end


"""
    bulk_generate_R(H, Ht, seedsize, Rsize, ntrial, genType; [dataset="HSBM", filefolder=homedir()*"/local-DHSG/data/", prob=prob])

Generate Rs for a hypergraph H. We randomly sample seedsize vertices as seeds,
then expand them to a R using method genType. We generate ntrial Rs.
"""
function bulk_generate_R(
    H::SparseMatrixCSC{Tf, Ti}, 
    Ht::SparseMatrixCSC{Tf, Ti},
    seedsize::Ti = 5,
    Rsize::Ti = 100,
    ntrial::Ti = 10,
    genType::String = "RW";
    dataset::String = "HSBM",
    filefolder = homedir()*"/local-DHSG/data/",
    prob::Tf = prob,
) where {Ti <: Integer, Tf <: AbstractFloat}
    Rs = Vector{Ti}[]     
    m, n = size(H)
    for i = 1:ntrial      
        seedset = sample(1:n, seedsize, replace=false)
        R = generate_R(H, Ht, seedset, Rsize, genType; prob=prob)
        push!(Rs, R)
    end
    savepath = filefolder * dataset * "/Rs-$seedsize-$Rsize-ntrial-$ntrial-$genType.mat"
    # save the results
    matwrite(savepath, Dict("Rs"=>Rs))
    return Rs
end


"""
    plot_runningtime(ϵs, results, labels, cols, savepath; [yscale=identity])

Plot the mean and stderr of the results.
"""
function plot_objective_ratio(
    ϵs::Vector{Tf},
    results::Vector{Any},
    labels,
    cols,
    savepath;
    yscale = identity,
) where{Tf}
    size_inches = (4, 3)
    size_pt = 90 .* size_inches
    fig = Figure(resolution = size_pt, fontsize = 16, figure_padding=1)
    axis = Axis(fig[1, 1], ylabel="ratio between objective values", yscale=yscale)
    for i = 1:length(labels)
        result = results[i] 
        label = labels[i]
        col = cols[i]
        result_mean = mean(result, dims=2)[:, 1]
        lines!(axis, ϵs, result_mean, label=label, color=(col, 1))
        scatter!(axis, ϵs, result_mean, color=(col, 1))
        n = size(result, 1)
        lowcurve = zeros(Tf, n)
        upcurve = zeros(Tf, n)
        for i = 1:n
            stderr = sem(result[i, :])
            low = result_mean[i] - stderr
            up = min(1.0, result_mean[i] + stderr)
            lowcurve[i] = low
            upcurve[i] = up
        end
        band!(axis, ϵs, upcurve, lowcurve, color = (col, 0.2))
    end
    axislegend(axis, position=:rb, merge = true, unique = true)
    save(savepath, fig)
end