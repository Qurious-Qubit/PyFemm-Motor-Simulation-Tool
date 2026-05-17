import os
import math
import csv
import numpy as np
from .fea_engine import run_base_fea_sweep

def execute_transient(
    *, motor_state, I_rms, gamma_elec_deg, theta_start_elec, theta_end_elec, num_steps, 
    Ns_rpm, active_filename, band_name, delete_files
):
    sim_folder = os.path.join(motor_state.base_folder, f"Transient_{I_rms}A_{gamma_elec_deg}deg_{Ns_rpm}RPM")
    if not os.path.exists(sim_folder): os.makedirs(sim_folder)
        
    print(f"\n--- Starting TRANSIENT Sweep ({num_steps} steps) on {active_filename} ---")
    
    res = run_base_fea_sweep(
        motor_state=motor_state, I_rms=I_rms, gamma_elec_deg=gamma_elec_deg, 
        theta_start_elec=theta_start_elec, theta_end_elec=theta_end_elec, num_steps=num_steps, 
        active_filename=active_filename, temp_folder=sim_folder, 
        band_name=band_name, delete_files=delete_files
    )
    
    # Dynamic Time Math (E_abc = dFlux/dt)
    omega_mech_rad_s = Ns_rpm * (2 * math.pi / 60.0)
    theta_mech_rad = np.radians(res["Theta_Mech_deg"])
    
    res["V_A_Volts"] = np.gradient(res["Flux_A_Wb"], theta_mech_rad) * omega_mech_rad_s
    res["V_B_Volts"] = np.gradient(res["Flux_B_Wb"], theta_mech_rad) * omega_mech_rad_s
    res["V_C_Volts"] = np.gradient(res["Flux_C_Wb"], theta_mech_rad) * omega_mech_rad_s
    
    # Dynamic DQ0 Voltages
    theta_rad = np.radians(res["Theta_Elec_deg"])
    v_a = np.array(res["V_A_Volts"]); v_b = np.array(res["V_B_Volts"]); v_c = np.array(res["V_C_Volts"])
    cos_t = np.cos(theta_rad); cos_t_120 = np.cos(theta_rad - 2*np.pi/3); cos_t_240 = np.cos(theta_rad + 2*np.pi/3)
    sin_t = np.sin(theta_rad); sin_t_120 = np.sin(theta_rad - 2*np.pi/3); sin_t_240 = np.sin(theta_rad + 2*np.pi/3)
    
    res["V_d_Volts"] = (2/3) * (v_a * cos_t + v_b * cos_t_120 + v_c * cos_t_240)
    res["V_q_Volts"] = -(2/3) * (v_a * sin_t + v_b * sin_t_120 + v_c * sin_t_240)
    
    # Save Files
    with open(os.path.join(sim_folder, "transient_waveforms.csv"), 'w', newline='') as f:
        w = csv.writer(f); w.writerow(res.keys()); w.writerows(zip(*[res[k] for k in res.keys()]))
        
    def rms(arr): return np.sqrt(np.mean(np.square(arr)))
    
    summary = {
        "I_RMS_A": rms(res["I_A"]), "V_RMS_Volts": rms(res["V_A_Volts"]),
        "Avg_Torque_Nm": np.mean(res["Torque_Nm"]), "Torque_Ripple_pk2pk": np.max(res["Torque_Nm"]) - np.min(res["Torque_Nm"]),
        "Avg_I_d_A": np.mean(res["I_d_A"]), "Avg_I_q_A": np.mean(res["I_q_A"]),
        "Avg_V_d_Volts": np.mean(res["V_d_Volts"]), "Avg_V_q_Volts": np.mean(res["V_q_Volts"]),
        "Avg_Flux_d_Wb": np.mean(res["Flux_d_Wb"]), "Avg_Flux_q_Wb": np.mean(res["Flux_q_Wb"])
    }
    with open(os.path.join(sim_folder, "transient_summary.csv"), 'w', newline='') as f:
        w = csv.writer(f); w.writerow(summary.keys()); w.writerow(summary.values())
        
    print(f"Transient Data Saved to: {sim_folder}")
    return res