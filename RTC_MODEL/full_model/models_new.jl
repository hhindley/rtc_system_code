function full_ODEs!(dydt, initial, params, t) # 3
    kr, L, c, sigma, atp, k1, k2, k3, k4, k5, w_rtc, theta_rtc, max, thr, d, k1_a, k2_a, k3_a, ka_b, kb_b, kc_b, k = params
    rm_a, rtca, rm_b, rtcb, rm_r, rtcr, rdrtca, rtrtcb, ribo_h, ribo_d, ribo_t = initial 
    drm_adt, drtcadt, drm_bdt, drtcbdt, drm_rdt, drtcrdt, drdrtcadt, drtrtcbdt, dribo_hdt, dribo_ddt, dribo_tdt = zeros(length(dydt))

    n = 6
    gr = 0*ribo_h 
    
    alpha = ribo_t/kr
    fa = (1+alpha)^n/(L*((1+c*alpha)^n)+(1+alpha)^n) # fraction of active RtcR
    ra = fa*rtcr # amount of active RtcR
    v = ra*k1*sigma*k3*atp*k5/(k1*sigma*k3*atp+k1*sigma*(k4+k5)+k2*(k4+k5)+k3*atp*k5) # rate of open complex formation
    
    # rtca
    tscr = v*(w_rtc*atp/(theta_rtc+atp)) # transcription rate 
    tlel = max*atp/(thr+atp) # translation elongation rate
    tlr = ribo_h*rm_a*tlel # translation rate so formation of Rtc 

    # rtcb
    tscr_b = v*(w_rtc*atp/(theta_rtc+atp)) # transcription rate 
    tlel_b = max*atp/(thr+atp) # translation elongation rate
    tlr_b = ribo_h*rm_b*tlel_b # translation rate so formation of Rtc 

    # rtcr 
    rtcr_tscr = w_rtc*atp/(theta_rtc+atp)
    rtcr_tlel = max*atp/(thr+atp)
    rtcr_tlr = ribo_h*rm_r*rtcr_tlel

    dydt[1] = tscr - gr*rm_a - d*rm_a
    dydt[2] = tlr - k1_a*rtca*ribo_d + k2_a*rdrtca + k3_a*atp*rdrtca - gr*rtca
    dydt[3] = tscr_b - gr*rm_b - d*rm_b
    dydt[4] = tlr_b - ka_b*rtcb*ribo_t + kb_b*rtrtcb + kc_b*atp*rtrtcb - gr*rtcb
    dydt[5] = rtcr_tscr - gr*rm_r - d*rm_r
    dydt[6] = rtcr_tlr - gr*rtcr
    dydt[7] = k1_a*rtca*ribo_d - k2_a*rdrtca - k3_a*atp*rdrtca
    dydt[8] = ka_b*rtcb*ribo_t - kb_b*rtrtcb - kc_b*atp*rtrtcb
    dydt[9] = kc_b*atp*rtrtcb - k*ribo_h - gr*ribo_h
    dydt[10] = k*ribo_h - k1_a*(rtca-rdrtca)*ribo_d + k2_a*rdrtca - gr*ribo_d
    dydt[11] = k3_a*atp*rdrtca - ka_b*(rtcb-rtrtcb)*ribo_t + kb_b*rtrtcb - gr*ribo_t

end

function mix1!(dydt, initial, params, t) # 1
    kr, L, c, sigma, atp, k1, k2, k3, k4, k5, w_rtc, theta_rtc, max, thr, d, k1_a, k2_a, k3_a, ka_b, kb_b, kc_b, k, rtca_tot, rtcb_tot = params
    rm_a, rtca, rm_b, rtcb, rm_r, rtcr, ribo_h, ribo_d, ribo_t = initial 
    drm_adt, drtcadt, drm_bdt, drtcbdt, drm_rdt, drtcrdt, dribo_hdt, dribo_ddt, dribo_tdt = zeros(length(dydt))

    n = 6


    # rtca_tot = rtca + rdrtca_0
    # rtcb_tot = rtcb + rtrtcb_0


    rtrtcb = ka_b*ribo_t*rtcb_tot/(ka_b*ribo_t+kb_b+kc_b*atp)
    rdrtca = k1_a*ribo_d*rtca_tot/(k1_a*ribo_d+k2_a+k3_a*atp)

    gr = 0*ribo_h 
    
    alpha = ribo_t/kr
    fa = (1+alpha)^n/(L*((1+c*alpha)^n)+(1+alpha)^n) # fraction of active RtcR
    ra = fa*rtcr # amount of active RtcR
    v = ra*k1*sigma*k3*atp*k5/(k1*sigma*k3*atp+k1*sigma*(k4+k5)+k2*(k4+k5)+k3*atp*k5) # rate of open complex formation
    
    # rtca
    tscr = v*(w_rtc*atp/(theta_rtc+atp)) # transcription rate 
    tlel = max*atp/(thr+atp) # translation elongation rate
    tlr = ribo_h*rm_a*tlel # translation rate so formation of Rtc 

    # rtcb
    tscr_b = v*(w_rtc*atp/(theta_rtc+atp)) # transcription rate 
    tlel_b = max*atp/(thr+atp) # translation elongation rate
    tlr_b = ribo_h*rm_b*tlel_b # translation rate so formation of Rtc 

    # rtcr 
    rtcr_tscr = w_rtc*atp/(theta_rtc+atp)
    rtcr_tlel = max*atp/(thr+atp)
    rtcr_tlr = ribo_h*rm_r*rtcr_tlel


    dydt[1] = tscr - gr*rm_a - d*rm_a
    dydt[2] = tlr - k1_a*rtca*ribo_d + k2_a*rdrtca + k3_a*atp*rdrtca - gr*rtca
    dydt[3] = tscr_b - gr*rm_b - d*rm_b
    dydt[4] = tlr_b - ka_b*rtcb*ribo_t + kb_b*rtrtcb + kc_b*atp*rtrtcb - gr*rtcb
    dydt[5] = rtcr_tscr - gr*rm_r - d*rm_r
    dydt[6] = rtcr_tlr - gr*rtcr

    dydt[7] = kc_b*atp*rtrtcb - k*ribo_h 
    dydt[8] = k*ribo_h - k1_a*rtca*ribo_d + k2_a*rdrtca 
    dydt[9] = k3_a*atp*rdrtca - ka_b*rtcb*ribo_t + kb_b*rtrtcb 
end

function mix2!(dydt, initial, params, t) # 4
    kr, L, c, sigma, atp, k1, k2, k3, k4, k5, w_rtc, theta_rtc, max, thr, d, k1_a, k2_a, k3_a, ka_b, kb_b, kc_b, k, ribo_tot = params
    rm_a, rtca, rm_b, rtcb, rm_r, rtcr, rdrtca, rtrtcb, ribo_d, ribo_t = initial 
    drm_adt, drtcadt, drm_bdt, drtcbdt, drm_rdt, drtcrdt, drdrtcadt, drtrtcbdt, dribo_ddt, dribo_tdt = zeros(length(dydt))

    n = 6

    ribo_h = ribo_tot - ribo_d - ribo_t - rdrtca - rtrtcb

    gr = 0*ribo_h 
    
    alpha = ribo_t/kr
    fa = (1+alpha)^n/(L*((1+c*alpha)^n)+(1+alpha)^n) # fraction of active RtcR
    ra = fa*rtcr # amount of active RtcR
    v = ra*k1*sigma*k3*atp*k5/(k1*sigma*k3*atp+k1*sigma*(k4+k5)+k2*(k4+k5)+k3*atp*k5) # rate of open complex formation
    
    # rtca
    tscr = v*(w_rtc*atp/(theta_rtc+atp)) # transcription rate 
    tlel = max*atp/(thr+atp) # translation elongation rate
    tlr = ribo_h*rm_a*tlel # translation rate so formation of Rtc 

    # rtcb
    tscr_b = v*(w_rtc*atp/(theta_rtc+atp)) # transcription rate 
    tlel_b = max*atp/(thr+atp) # translation elongation rate
    tlr_b = ribo_h*rm_b*tlel_b # translation rate so formation of Rtc 

    # rtcr 
    rtcr_tscr = w_rtc*atp/(theta_rtc+atp)
    rtcr_tlel = max*atp/(thr+atp)
    rtcr_tlr = ribo_h*rm_r*rtcr_tlel

    dydt[1] = tscr - gr*rm_a - d*rm_a
    dydt[2] = tlr - k1_a*rtca*ribo_d + k2_a*rdrtca + k3_a*atp*rdrtca - gr*rtca
    dydt[3] = tscr_b - gr*rm_b - d*rm_b
    dydt[4] = tlr_b - ka_b*rtcb*ribo_t + kb_b*rtrtcb + kc_b*atp*rtrtcb - gr*rtcb
    dydt[5] = rtcr_tscr - gr*rm_r - d*rm_r
    dydt[6] = rtcr_tlr - gr*rtcr

    dydt[7] = k1_a*rtca*ribo_d - k2_a*rdrtca - k3_a*atp*rdrtca
    dydt[8] = ka_b*rtcb*ribo_t - kb_b*rtrtcb - kc_b*atp*rtrtcb

    dydt[9] = k*ribo_h - k1_a*rtca*ribo_d + k2_a*rdrtca
    dydt[10] = k3_a*atp*rdrtca - ka_b*rtcb*ribo_t + kb_b*rtrtcb
end

function full_alg!(dydt, initial, params, t) # 2
    kr, L, c, sigma, atp, k1, k2, k3, k4, k5, w_rtc, theta_rtc, max, thr, d, k1_a, k2_a, k3_a, ka_b, kb_b, kc_b, k, rtca_tot, rtcb_tot, ribo_tot = params
    rm_a, rtca, rm_b, rtcb, rm_r, rtcr, ribo_d, ribo_t = initial 
    drm_adt, drtcadt, drm_bdt, drtcbdt, drm_rdt, drtcrdt, dribo_ddt, dribo_tdt = zeros(length(dydt))

    n = 6


    # rtca_tot = rtca + rdrtca_0
    # rtcb_tot = rtcb + rtrtcb_0


    rtrtcb = ka_b*ribo_t*rtcb_tot/(ka_b*ribo_t+kb_b+kc_b*atp)

    rdrtca = k1_a*ribo_d*rtca_tot/(k1_a*ribo_d+k2_a+k3_a*atp)

    ribo_h = ribo_tot - ribo_d - ribo_t - rdrtca - rtrtcb

    gr = 0*ribo_h
    
    alpha = ribo_t/kr
    fa = (1+alpha)^n/(L*((1+c*alpha)^n)+(1+alpha)^n) # fraction of active RtcR
    ra = fa*rtcr # amount of active RtcR
    v = ra*k1*sigma*k3*atp*k5/(k1*sigma*k3*atp+k1*sigma*(k4+k5)+k2*(k4+k5)+k3*atp*k5) # rate of open complex formation
    
    # rtca
    tscr = v*(w_rtc*atp/(theta_rtc+atp)) # transcription rate 
    tlel = max*atp/(thr+atp) # translation elongation rate
    tlr = ribo_h*rm_a*tlel # translation rate so formation of Rtc 

    # rtcb
    tscr_b = v*(w_rtc*atp/(theta_rtc+atp)) # transcription rate 
    tlel_b = max*atp/(thr+atp) # translation elongation rate
    tlr_b = ribo_h*rm_b*tlel_b # translation rate so formation of Rtc 

    # rtcr 
    rtcr_tscr = w_rtc*atp/(theta_rtc+atp)
    rtcr_tlel = max*atp/(thr+atp)
    rtcr_tlr = ribo_h*rm_r*rtcr_tlel

    dydt[1] = tscr - gr*rm_a - d*rm_a
    dydt[2] = tlr - k1_a*rtca*ribo_d + k2_a*rdrtca + k3_a*atp*rdrtca - gr*rtca
    dydt[3] = tscr_b - gr*rm_b - d*rm_b
    dydt[4] = tlr_b - ka_b*rtcb*ribo_t + kb_b*rtrtcb + kc_b*atp*rtrtcb - gr*rtcb
    dydt[5] = rtcr_tscr - gr*rm_r - d*rm_r
    dydt[6] = rtcr_tlr - gr*rtcr

    dydt[7] = k*ribo_h - k1_a*rtca*ribo_d + k2_a*rdrtca
    dydt[8] = k3_a*atp*rdrtca - ka_b*rtcb*ribo_t + kb_b*rtrtcb
end



function RTC_SYSTEM!(dydt, initial, params, t) # 2
    kr, L, c, sigma, atp, k1, k2, k3, k4, k5, w_rtc, theta_rtc, max, thr, d, k1_a, k2_a, k3_a, ka_b, kb_b, kc_b, k, ribo_tot = params
    rm_a, rtca, rm_b, rtcb, rm_r, rtcr, ribo_d, ribo_t = initial 
    drm_adt, drtcadt, drm_bdt, drtcbdt, drm_rdt, drtcrdt, dribo_ddt, dribo_tdt = zeros(length(dydt))

    n = 6
    
    rdrtca = k1_a*ribo_d*rtca/(k1_a*ribo_d+k2_a+k3_a*atp)
    rtrtcb = ka_b*ribo_t*rtcb/(ka_b*ribo_t+kb_b+kc_b*atp)
    ribo_h = ribo_tot - ribo_d - ribo_t - rdrtca - rtrtcb

    gr = 0*ribo_h
    
    alpha = ribo_t/kr
    fa = (1+alpha)^n/(L*((1+c*alpha)^n)+(1+alpha)^n) # fraction of active RtcR
    ra = fa*rtcr # amount of active RtcR
    v = ra*k1*sigma*k3*atp*k5/(k1*sigma*k3*atp+k1*sigma*(k4+k5)+k2*(k4+k5)+k3*atp*k5) # rate of open complex formation
    
    # rtca
    tscr = v*(w_rtc*atp/(theta_rtc+atp)) # transcription rate 
    tlel = max*atp/(thr+atp) # translation elongation rate
    tlr = ribo_h*rm_a*tlel # translation rate so formation of Rtc 

    # rtcb
    tscr_b = v*(w_rtc*atp/(theta_rtc+atp)) # transcription rate 
    tlel_b = max*atp/(thr+atp) # translation elongation rate
    tlr_b = ribo_h*rm_b*tlel_b # translation rate so formation of Rtc 

    # rtcr 
    rtcr_tscr = w_rtc*atp/(theta_rtc+atp)
    rtcr_tlel = max*atp/(thr+atp)
    rtcr_tlr = ribo_h*rm_r*rtcr_tlel

    dydt[1] = tscr - gr*rm_a - d*rm_a
    dydt[2] = tlr - gr*rtca
    dydt[3] = tscr_b - gr*rm_b - d*rm_b
    dydt[4] = tlr_b - gr*rtcb
    dydt[5] = rtcr_tscr - gr*rm_r - d*rm_r
    dydt[6] = rtcr_tlr - gr*rtcr

    dydt[7] = k*ribo_h - k3_a*rdrtca - gr*ribo_d
    dydt[8] = k3_a*rdrtca - kc_b*rtrtcb - gr*ribo_t
    
end

function new_ODEs!(dydt, initial, params, t) # 3
    kr, L, c, sigma, atp, k1, k2, k3, k4, k5, w_rtc, theta_rtc, max, thr, d, k1_a, k2_a, k3_a, ka_b, kb_b, kc_b, k = params
    rm_a, rtca, rm_b, rtcb, rm_r, rtcr, rdrtca, rtrtcb, ribo_h, ribo_d, ribo_t = initial 
    drm_adt, drtcadt, drm_bdt, drtcbdt, drm_rdt, drtcrdt, drdrtcadt, drtrtcbdt, dribo_hdt, dribo_ddt, dribo_tdt = zeros(length(dydt))

    n = 6
    gr = 0*ribo_h 
    
    alpha = ribo_t/kr
    fa = (1+alpha)^n/(L*((1+c*alpha)^n)+(1+alpha)^n) # fraction of active RtcR
    ra = fa*rtcr # amount of active RtcR
    v = ra*k1*sigma*k3*atp*k5/(k1*sigma*k3*atp+k1*sigma*(k4+k5)+k2*(k4+k5)+k3*atp*k5) # rate of open complex formation
    
    # rtca
    tscr = v*(w_rtc*atp/(theta_rtc+atp)) # transcription rate 
    tlel = max*atp/(thr+atp) # translation elongation rate
    tlr = ribo_h*rm_a*tlel # translation rate so formation of Rtc 

    # rtcb
    tscr_b = v*(w_rtc*atp/(theta_rtc+atp)) # transcription rate 
    tlel_b = max*atp/(thr+atp) # translation elongation rate
    tlr_b = ribo_h*rm_b*tlel_b # translation rate so formation of Rtc 

    # rtcr 
    rtcr_tscr = w_rtc*atp/(theta_rtc+atp)
    rtcr_tlel = max*atp/(thr+atp)
    rtcr_tlr = ribo_h*rm_r*rtcr_tlel

    dydt[1] = tscr - gr*rm_a - d*rm_a
    dydt[2] = tlr - gr*rtca # rtca_tot
    dydt[3] = tscr_b - gr*rm_b - d*rm_b
    dydt[4] = tlr_b - gr*rtcb # rtcb_tot
    dydt[5] = rtcr_tscr - gr*rm_r - d*rm_r
    dydt[6] = rtcr_tlr - gr*rtcr

    # rtca_tot = rtca+rdrtca
    dydt[7] = k1_a*(rtca-rdrtca)*ribo_d - k2_a*rdrtca - k3_a*atp*rdrtca 
    dydt[8] = ka_b*(rtcb-rtrtcb)*ribo_t - kb_b*rtrtcb - kc_b*atp*rtrtcb

    dydt[9] = kc_b*rtrtcb - k*ribo_h - gr*ribo_h
    dydt[10] = k*ribo_h - k3_a*rdrtca - gr*ribo_d
    dydt[11] = k3_a*rdrtca - kc_b*rtrtcb - gr*ribo_t

end