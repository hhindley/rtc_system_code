using Plots, DifferentialEquations, DataFrames, Statistics, Measures
include("fullmodel_functions.jl")
include("fullmodel_params.jl")
include("models_new.jl")
# plotlyjs()

sol = simple_solve!(mix2!, init_mix2, tspan, params_mix2)
plot(sol[2:end], xaxis=(:log10, (1,Inf)), yaxis=(:log10, (1,Inf)), labels=labels_mix2, palette=:tab10)


# splitting species and looking at solution results
rm_a, rtca, rm_b, rtcb, rm_r, rtcr, ribo_h, ribo_d, ribo_t, rdrtca, rtrtcb = split_mix2!(sol, species_mix2)
plot(sol.t, [rm_a, rtca, rm_b, rtcb, rm_r, rtcr, rdrtca, rtrtcb, ribo_d, ribo_h, ribo_t], xaxis=(:log10, (1,Inf)), yaxis=(:log10, (1,Inf)), labels=labels_all)

p_sol = plot(sol.t, [ribo_d+ribo_h+ribo_t+rdrtca+rtrtcb], xaxis=(:log10, (1,Inf)), margin=16mm)

p_rtc = plot(sol.t, xlabel="t", [rtca, rtcb, rtcr], labels=["RtcA" "RtcB" "RtcR"], xaxis=(:log10, (1,Inf)), yaxis=(:log10, (1,Inf)))
p_mrna = plot(sol.t, xlabel="t", [rm_a, rm_b, rm_r], labels=["RtcA-mRNA" "RtcB-mRNA" "RtcR-mRNA"], xaxis=(:log10, (1,Inf)))
p_ribo = plot(sol.t, xlabel="t", [ribo_d, ribo_h, ribo_t], labels=["Rd" "Rh" "Rt"], xaxis=(:log10, (1,Inf)))#, title="Ribosomes")
p_int = plot(sol.t, xlabel="t", [rdrtca, rtrtcb], title="Intermediates", labels=["RdRtcA" "RtRtcB"], xaxis=(:log10, (1,Inf)))

plot(p_sol, p_rtc, p_mrna, p_ribo, p_int, size=(1000,1000))

plot(sol.t, ribo_d, xlabel="t", title="Damaged ribosomes", legend=false, xaxis=(:log10, (1,Inf)))
plot(sol.t, ribo_t, xlabel="t", title="Tagged ribosomes", legend=false, xaxis=(:log10, (1,Inf)))
plot(sol.t, ribo_h, xlabel="t", title="Healthy ribosomes", legend=false, xaxis=(:log10, (1,Inf)))


# varying parameters
rd_range = range(0, 10, length=6)
include("fullmodel_params.jl")
rd_p1, rd_p2 = logplots(rd_range, init_mix2, 9, mix2!, init_mix2, params_mix2, species_mix2, labels_mix2, "Rd")
rd_p1
rd_p2

atp_range = range(0, 20, length=6) 
include("fullmodel_params.jl")
sols, sols_t, rm_a_sols, rtca_sols, rm_b_sols, rtcb_sols, rm_r_sols, rtcr_sols, ribo_h_sols, ribo_d_sols, ribo_t_sols, rdrtca_sols, rtrtcb_sols, ribo_tot_sols = vary_param(atp_range, params_rtc, 5, RTC_SYSTEM!, init_rtc, params_rtc, species_rtc)
plot(sols_t[3], [rm_a_sols[3], rtca_sols[3], rm_b_sols[3], rtcb_sols[3], rm_r_sols[3], rtcr_sols[3], rdrtca_sols[3], rtrtcb_sols[3], ribo_d_sols[3], ribo_h_sols[3], ribo_t_sols[3]], xlabel="t", yaxis=(:log10, (1,Inf)), xaxis=(:log10, (1,Inf)), labels=labels_all, palette=:tab10, title="5 - new model")
plot(sols_t[3], [rdrtca_sols[3], rtrtcb_sols[3]], xaxis=(:log10, (1,Inf)), labels=["RdRtcA" "RtRtcB"])


atp_p1, atp_p2 = logplots(atp_range, params_mix2, 5, mix2!, init_mix2, params_mix2, species_mix2, labels_mix2, "ATP")
atp_p1
atp_p2

kr_range = range(1, 6, length=6)
include("fullmodel_params.jl") 
kr_p1, kr_p2 = logplots(kr_range, params_mix2, 1, mix2!, init_mix2, params_mix2, species_mix2, labels_mix2, "kr")
kr_p1
kr_p2

L_range = range(0, 10, length=6) 
include("fullmodel_params.jl")
L_p1, L_p2 = logplots(L_range, params_mix2, 2, mix2!, init_mix2, params_mix2, species_mix2, labels_mix2, "L")
L_p1
L_p2

c_range = range(0, 100, length=6)
include("fullmodel_params.jl")
c_p1, c_p2 = logplots(c_range, params_mix2, 3, mix2!, init_mix2, params_mix2, species_mix2, labels_mix2, "c")
c_p1
c_p2

sigma_range = range(0, 10, length=6)
include("fullmodel_params.jl")
sigma_p1, sigma_p2 = logplots(sigma_range, params_mix2, 4, mix2!, init_mix2, params_mix2, species_mix2, labels_mix2, "sigma")
sigma_p1
sigma_p2

k1_range = range(0, 2, length=6)
include("fullmodel_params.jl")
k1_p1, k1_p2 = logplots(k1_range, params_mix2, 6, mix2!, init_mix2, params_mix2, species_mix2, labels_mix2, "k1")
k1_p1
k1_p2

k2_range = range(0, 2, length=6)
include("fullmodel_params.jl")
k2_p1, k2_p2 = logplots(k2_range, params_mix2, 7, mix2!, init_mix2, params_mix2, species_mix2, labels_mix2, "k2")
k2_p1
k2_p2

k3_range = range(0, 2, length=6)
include("fullmodel_params.jl")
k3_p1, k3_p2 = logplots(k3_range, params_mix2, 8, mix2!, init_mix2, params_mix2, species_mix2, labels_mix2, "k3")
k3_p1
k3_p2

k4_range = range(0, 2, length=6)
include("fullmodel_params.jl")
k4_p1, k4_p2 = logplots(k4_range, params_mix2, 9, mix2!, init_mix2, params_mix2, species_mix2, labels_mix2, "k4")
k4_p1
k4_p2

k5_range = range(0, 2, length=6)
include("fullmodel_params.jl")
k5_p1, k5_p2 = logplots(k5_range, params_mix2, 10, mix2!, init_mix2, params_mix2, species_mix2, labels_mix2, "k5")
k5_p1
k5_p2

wrtc_range = range(0, 10, length=6)
include("fullmodel_params.jl")
wrtc_p1, wrtc_p2 = logplots(wrtc_range, params_mix2, 11, mix2!, init_mix2, params_mix2, species_mix2, labels_mix2, "wrtc")
wrtc_p1
wrtc_p2

thrtc_range = range(0, 20, length=6)
include("fullmodel_params.jl")
thrtc_p1, thrtc_p2 = logplots(thrtc_range, params_mix2, 12, mix2!, init_mix2, params_mix2, species_mix2, labels_mix2, "thrtc")
thrtc_p1
thrtc_p2

max_range = range(0, 10, length=6)
include("fullmodel_params.jl")
max_p1, max_p2 = logplots(max_range, params_mix2, 13, mix2!, init_mix2, params_mix2, species_mix2, labels_mix2, "max")
max_p1
max_p2

thr_range = range(0, 10, length=6)
include("fullmodel_params.jl")
thr_p1, thr_p2 = logplots(thr_range, params_mix2, 14, mix2!, init_mix2, params_mix2, species_mix2, labels_mix2, "thr")
thr_p1
thr_p2

d_range = range(0, 5, length=6)
include("fullmodel_params.jl")
d_p1, d_p2 = logplots(d_range, params_mix2, 15, mix2!, init_mix2, params_mix2, species_mix2, labels_mix2, "d")
d_p1
d_p2

k1a_range = range(0, 2, length=6)
include("fullmodel_params.jl")
k1a_p1, k1a_p2 = logplots(k1a_range, params_mix2, 16, mix2!, init_mix2, params_mix2, species_mix2, labels_mix2, "k1a")
k1a_p1
k1a_p2

k2a_range = range(0, 2, length=6)
include("fullmodel_params.jl")
k2a_p1, k2a_p2 = logplots(k2a_range, params_mix2, 17, mix2!, init_mix2, params_mix2, species_mix2, labels_mix2, "k2a")
k2a_p1
k2a_p2

k3a_range = range(0, 2, length=6)
include("fullmodel_params.jl")
k3a_p1, k3a_p2 = logplots(k3a_range, params_mix2, 18, mix2!, init_mix2, params_mix2, species_mix2, labels_mix2, "k3a")
k3a_p1
k3a_p2

kab_range = range(0, 2, length=6)
include("fullmodel_params.jl")
kab_p1, kab_p2 = logplots(kab_range, params_mix2, 19, mix2!, init_mix2, params_mix2, species_mix2, labels_mix2, "kab")
kab_p1
kab_p2

kbb_range = range(0, 2, length=6)
include("fullmodel_params.jl")
kbb_p1, kbb_p2 = logplots(kbb_range, params_mix2, 20, mix2!, init_mix2, params_mix2, species_mix2, labels_mix2, "kbb")
kbb_p1
kbb_p2

kcb_range = range(0, 2, length=6)
include("fullmodel_params.jl")
kcb_p1, kcb_p2 = logplots(kcb_range, params_mix2, 21, mix2!, init_mix2, params_mix2, species_mix2, labels_mix2, "kcb")
kcb_p1
kcb_p2

k_range = range(0, 2, length=6)
include("fullmodel_params.jl")
k_p1, k_p2 = logplots(k_range, params_mix2, 22, mix2!, init_mix2, params_mix2, species_mix2, labels_mix2, "k")
k_p1
k_p2

rtca_range = range(0, 20, length=6) # over 10 and start seeing negative values?
include("fullmodel_params.jl")
rtca_p1, rtca_p2 = logplots(rtca_range, init_mix2, 2, mix2!, init_mix2, params_mix2, species_mix2, labels_mix2, "rtca")
rtca_p1
rtca_p2

rtcb_range = range(0, 20, length=6) # over 10 and start seeing negative values?
include("fullmodel_params.jl")
rtcb_p1, rtcb_p2 = logplots(rtcb_range, init_mix2, 4, mix2!, init_mix2, params_mix2, species_mix2, labels_mix2, "rtcb")
rtcb_p1
rtcb_p2

rt_range = range(0, 10, length=6) # over 10 and start seeing negative values - because of the equation for rt - change rt in full model to ode and see if that works (probs will) 
include("fullmodel_params.jl")
rt_p1, rt_p2 = logplots(rt_range, init_mix2, 10, mix2!, init_mix2, params_mix2, species_mix2, labels_mix2, "Rt")
rt_p1
rt_p2






