using CSV
using DataFrames
using Plots
using Distributions

Tgvalues = [25.4,45.1,102.1]                                                    # Tg values
wpInitsArray = [0.783,0.822,0.858,0.887,0.893,0.909,0.919,0.926,0.936]          # wp initialization for each layer

for (ntg,tgv) in enumerate(Tgvalues)                                            # Iterate over TgValues
    v_read = CSV.read("data/paper/"*string(Tgvalues[ntg])*".csv", DataFrame; header=false) |> Matrix
    plt = Plots.plot()
    for (nwp,wpInit) in enumerate(wpInitsArray)                                     # Iterate over wp : layer increase with wp
        v = v_read[nwp,:]
        x = range(0, 83, length=10000)
        Plots.plot!(plt,x,pdf.(fit(Gamma, v),x), label = wpInit)
    end
    savefig(plt,"plots/MlayerChgWp/Paper/DeepestHistogram/alltg-"*string(tgv)*".png")
end
