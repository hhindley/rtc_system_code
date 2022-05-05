using Plots, DifferentialEquations, DataFrames, Statistics, Measures
include("full_model/fullmodel_functions.jl")
include("full_model/fullmodel_params.jl")
include("full_model/models_new.jl")

prob = ODEProblem(new_ODEs!, init_new, tspan, params_new)
sol = solve(prob, Rodas4())
p = plot(sol[2:end], yaxis=(:log10, (1,Inf)), xaxis=(:log10, (1,Inf)), labels=labels_new, title="full ODEs");

prob1 = ODEProblem(RTC_SYSTEM!, init_rtc, tspan, params_rtc)
sol1 = solve(prob1, Rodas4())
p1 = plot(sol1[2:end], yaxis=(:log10, (1,Inf)), xaxis=(:log10, (1,Inf)), labels=labels_rtc, title="reduced model");

plot(p, p1, size=(800,500))

rm_a, rtca, rm_b, rtcb, rm_r, rtcr, rdrtca, rtrtcb, ribo_d, ribo_h, ribo_t = split_full_model(species_full, sol)
rm_a1, rtca1, rm_b1, rtcb1, rm_r1, rtcr1, rdrtca1, rtrtcb1, ribo_d1, ribo_h1, ribo_t1 = split_reduced_model(species_rtc, sol1)

p_full = plot(sol.t, [rm_a, rtca, rm_b, rtcb, rm_r, rtcr, rdrtca, rtrtcb, ribo_d, ribo_h, ribo_t], yaxis=(:log10, (1,Inf)), xaxis=(:log10, (1,Inf)), labels=labels_new, title="full ODEs");
p1_full = plot(sol1.t, [rm_a1, rtca1, rm_b1, rtcb1, rm_r1, rtcr1, rdrtca1, rtrtcb1, ribo_d1, ribo_h1, ribo_t1], yaxis=(:log10, (1,Inf)), xaxis=(:log10, (1,Inf)), labels=labels_new, title="reduced model");
plot(p_full, p1_full, size=(800,500))

plot(sol.t, [rm_a, rm_b, rm_r], yaxis=(:log10, (1,Inf)), xaxis=(:log10, (1,Inf)), labels=["mRNA-A-full" "mRNA-B-full" "mRNA-R-full"])
plot!(sol1.t, [rm_a1, rm_b1, rm_r1], yaxis=(:log10, (1,Inf)), xaxis=(:log10, (1,Inf)), labels=["mRNA-A" "mRNA-B" "mRNA-R"])

plot(sol.t, [rtca, rtcb, rtcr], yaxis=(:log10, (1,Inf)), xaxis=(:log10, (1,Inf)), labels=["RtcA-full" "RtcB-full" "RtcR-full"])
plot!(sol1.t, [rtca1, rtcb1, rtcr1], yaxis=(:log10, (1,Inf)), xaxis=(:log10, (1,Inf)), labels=["RtcA" "RtcB" "RtcR"])

plot(sol.t, [rdrtca, rtrtcb], xaxis=(:log10, (1,Inf)), labels=["RdRtcA-full" "RtRtcB-full"])
plot!(sol1.t, [rdrtca1, rtrtcb1], xaxis=(:log10, (1,Inf)), labels=["RdRtcA" "RtRtcB"])

plot(sol.t, ribo_d, xaxis=(:log10, (1,Inf)), labels="Rd-full")
plot!(sol1.t, ribo_d1, xaxis=(:log10, (1,Inf)), labels="Rd")

plot(sol.t, ribo_h, xaxis=(:log10, (1,Inf)), labels="Rh-full")
plot!(sol1.t, ribo_h1, xaxis=(:log10, (1,Inf)), labels="Rh")

plot(sol.t, ribo_t, xaxis=(:log10, (1,Inf)), labels="Rt-full")
plot!(sol1.t, ribo_t1, xaxis=(:log10, (1,Inf)), labels="Rt")

plot(sol.t, [ribo_d, ribo_h, ribo_t], xaxis=(:log10, (1,Inf)), labels=["Rd-full" "Rh-full" "Rt-full"])
plot!(sol1.t, [ribo_d1, ribo_h1, ribo_t1], xaxis=(:log10, (1,Inf)), labels=["Rd" "Rh" "Rt"])

plot(sol.t, ribo_d+ribo_h+ribo_t+rdrtca+rtrtcb, xaxis=(:log10, (1,Inf)), labels="Rtot-full")
plot(sol1.t, ribo_d1+ribo_h1+ribo_t1+rdrtca1+rtrtcb1, xaxis=(:log10, (1,Inf)), margin=16mm, labels="Rtot")
