full_name_pref = "g_report_ekrt_like"
mid_name_pref = "g_report_midrap_ekrt_like"
K_factors = [2, 4]
M_factors = [0.1, 0.2, 0.3, 0.4, 0.5, 0.65, 0.8, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 20]
print("{} {} {} {} {} {} {} {} {} {}".format("K", "M", 
      "midrap_sat_p", "f_sat_p", "midrap_sat_n", "f_sat_n", 
      "midrap_sat_mc_p", "f_sat_mc_p", "midrap_sat_mc_n", "f_sat_mc_n"
    ))
for K in K_factors:
  for M in M_factors:
    full_name = full_name_pref+"_K="+str(K)+"_M="+str(M)
    mid_name = mid_name_pref+"_K="+str(K)+"_M="+str(M)

    with open(mid_name+"_SAT.dat", "r") as f:
      m_s_p = f.read().split('\n')[88].split(' ')[-1]
    with open(full_name+"_SAT.dat", "r") as f:
      f_s_p = f.read().split('\n')[88].split(' ')[-1]
    with open(mid_name+"_nPDF_SAT.dat", "r") as f:
      m_s_n = f.read().split('\n')[88].split(' ')[-1]
    with open(full_name+"_nPDF_SAT.dat", "r") as f:
      f_s_n = f.read().split('\n')[88].split(' ')[-1]
    with open(mid_name+"_SAT_MC.dat", "r") as f:
      m_s_mc_p = f.read().split('\n')[88].split(' ')[-1]
    with open(full_name+"_SAT_MC.dat", "r") as f:
      f_s_mc_p = f.read().split('\n')[88].split(' ')[-1]
    with open(mid_name+"_nPDF_SAT_MC.dat", "r") as f:
      m_s_mc_n = f.read().split('\n')[88].split(' ')[-1]
    with open(full_name+"_nPDF_SAT_MC.dat", "r") as f:
      f_s_mc_n = f.read().split('\n')[88].split(' ')[-1]
      
    
    print("{} {} {} {} {} {} {} {} {} {}".format(K, M, 
      m_s_p, f_s_p, m_s_n, f_s_n, 
      m_s_mc_p, f_s_mc_p, m_s_mc_n, f_s_mc_n
    ))