include("../header.jl")
using GMT
using HTTP

regions_gmt = gmtread("/vsicurl/https://github.com/datasets/geo-countries/raw/master/data/countries.geojson");
regions = CSV.File(homedir()*"/local-DHSG/src/Exp3-web-edu/all.csv") |> DataFrame;


function iso_a2_to_a3(code, regions)
      for i = 1:size(regions, 1)
           if code == regions[i, "alpha-2"]         
                  return regions[i, "alpha-3"] 
           end
      end
      @assert false @show code
end 


function dict2df(a, regions)
      region_codes = []
      counts = []
      for (key, val) in a
            key = uppercase(key)
            key = iso_a2_to_a3(key, regions)
            if val == 0
                  continue
            end
            push!(counts, 1.0*val)
            push!(region_codes, uppercase(key))
      end
      return DataFrame(Region=region_codes, Count=counts)
end


function normalized(a, b)
      c = deepcopy(a)
      for (key, val) in a
            c[key] = a[key] / b[key] 
      end
      return c 
end

res = matread(homedir()*"/local-DHSG/results/webgraph/webgraph-edu-ac.mat")
cn_region_counts = res["cn_region_counts"]
uk_region_counts = res["uk_region_counts"]
inter_region_counts = res["inter_region_counts"]
total_region_counts = res["total_region_counts"]


cn_region_counts_norm = normalized(cn_region_counts, total_region_counts)
uk_region_counts_norm = normalized(uk_region_counts, total_region_counts)
inter_region_counts_norm = normalized(inter_region_counts, total_region_counts)


cn_region_counts_norm["ca"] = cn_region_counts_norm["us"]
uk_region_counts_norm["ca"] = uk_region_counts_norm["us"]
inter_region_counts_norm["ca"] = inter_region_counts_norm["us"]


cn_df = dict2df(cn_region_counts_norm, regions)
uk_df = dict2df(uk_region_counts_norm, regions)
it_df = dict2df(inter_region_counts_norm, regions)
cn_zvals = polygonlevels(regions_gmt, string.(cn_df[!,1]), Float64.(cn_df[!,2]), att="ISO_A3");
uk_zvals = polygonlevels(regions_gmt, string.(uk_df[!,1]), Float64.(uk_df[!,2]), att="ISO_A3");
it_zvals = polygonlevels(regions_gmt, string.(it_df[!,1]), Float64.(it_df[!,2]), att="ISO_A3");

C = makecpt(range=(0.0, 0.1), extend=true, bg=true); 
figpath = homedir()*"/local-DHSG/figs/webgraph/"

GMT.plot(regions_gmt, region=(-180,180,-60,85), level=cn_zvals, cmap=C, proj=:guess, frame=:none, pen=(0.1, 200))
colorbar!(pos=(anchor=:LB,length=(4,0.3), offset=(-0.5,-5)), color=C,
     axes=(annot=0.05,), savefig=figpath*"cn.pdf")

GMT.plot(regions_gmt, region=(-180,180,-60,85), level=uk_zvals, cmap=C, proj=:guess, frame=:none, pen=(0.1, 200))
colorbar!(pos=(anchor=:LB,length=(4,0.3), offset=(-0.5,-5)), color=C,
     axes=(annot=0.05,), savefig=figpath*"uk.pdf")

GMT.plot(regions_gmt, region=(-180,180,-60,85), level=it_zvals, cmap=C, proj=:guess, frame=:none, pen=(0.1, 200))
colorbar!(pos=(anchor=:LB,length=(4,0.3), offset=(-0.5,-5)), color=C,
     axes=(annot=0.05,), savefig=figpath*"it.pdf")



            




