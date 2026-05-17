import os
import csv
import numpy as np
import math
from .fea_engine import run_base_fea_sweep

def execute_mtpa(
    *, motor_state, I_rms_list, theta_start_elec, theta_end_elec, num_steps, 
    Ns_rpm, coarse_step_deg, run_fine_search, fine_window_deg, fine_step_deg,
    active_filename, band_name, delete_files
):
    mtpa_root = os.path.join(motor_state.base_folder, "MTPA_Analysis")
    if not os.path.exists(mtpa_root): os.makedirs(mtpa_root)
    
    print(f"\n--- Starting MTPA Sweep on {active_filename} ---")
    
    omega_e = (Ns_rpm * 2 * math.pi / 60.0) * (motor_state.p / 2)
    
    master_mtpa_data = {"I_rms_A": [], "MTPA_Gamma_deg": [], "Avg_Torque_Nm": [], "Torque_Ripple_pk2pk": [], "I_d_A": [], "I_q_A": [], "Flux_d_Wb": [], "Flux_q_Wb": [], "V_d_Volts": [], "V_q_Volts": []}
    grid_map_data = {"I_rms_A": [], "Gamma_deg": [], "Avg_Torque_Nm": [], "Torque_Ripple_pk2pk": [], "I_d_A": [], "I_q_A": [], "Flux_d_Wb": [], "Flux_q_Wb": [], "V_d_Volts": [], "V_q_Volts": []}
    all_raw_data = [] 

    for current in I_rms_list:
        print(f"\n{'='*40}\nMTPA SEARCH: I_RMS = {current} A\n{'='*40}")
        current_folder = os.path.join(mtpa_root, f"I_{current}A")
        
        # --- Generate Gated Coarse Gammas ---
        coarse_gammas = list(np.arange(0, 90 + coarse_step_deg, coarse_step_deg))
        if 90 not in coarse_gammas: coarse_gammas.append(90)
        coarse_gammas = sorted(list(set([g for g in coarse_gammas if g <= 90])))
        
        def evaluate_gamma(g):
            g_folder = os.path.join(current_folder, f"Gamma_{g:02g}deg")
            if not os.path.exists(g_folder): os.makedirs(g_folder)
            
            res = run_base_fea_sweep(
                motor_state=motor_state, I_rms=current, gamma_elec_deg=g, 
                theta_start_elec=theta_start_elec, theta_end_elec=theta_end_elec, num_steps=num_steps, 
                active_filename=active_filename, temp_folder=g_folder, 
                band_name=band_name, delete_files=delete_files
            )
            
            with open(os.path.join(g_folder, "simulation_results.csv"), 'w', newline='') as f:
                w = csv.writer(f); w.writerow(res.keys()); w.writerows(zip(*[res[k] for k in res.keys()]))
            
            avg_tq = np.mean(res["Torque_Nm"][:-1]); ripple = np.max(res["Torque_Nm"]) - np.min(res["Torque_Nm"])
            avg_id = np.mean(res["I_d_A"][:-1]); avg_iq = np.mean(res["I_q_A"][:-1])
            avg_fd = np.mean(res["Flux_d_Wb"][:-1]); avg_fq = np.mean(res["Flux_q_Wb"][:-1])
            v_d = -omega_e * avg_fq; v_q = omega_e * avg_fd
            
            for k, v in zip(grid_map_data.keys(), [current, g, avg_tq, ripple, avg_id, avg_iq, avg_fd, avg_fq, v_d, v_q]):
                grid_map_data[k].append(v)
            all_raw_data.append({"current": current, "gamma": g, "res": res})
            
            return avg_tq, [current, g, avg_tq, ripple, avg_id, avg_iq, avg_fd, avg_fq, v_d, v_q]

        # 1. Coarse Sweep execution
        best_coarse_torque = -float('inf'); best_coarse_gamma = 0
        for g in coarse_gammas:
            print(f" -> Coarse: Gamma = {g}°")
            tq, _ = evaluate_gamma(g)
            if tq > best_coarse_torque:
                best_coarse_torque = tq; best_coarse_gamma = g

        best_res_array = None

        # 2. Fine Sweep execution (If enabled)
        if run_fine_search:
            min_g = max(0, best_coarse_gamma - (fine_window_deg / 2))
            max_g = min(90, best_coarse_gamma + (fine_window_deg / 2))
            
            fine_gammas = list(np.arange(min_g, max_g + fine_step_deg, fine_step_deg))
            if max_g == 90 and 90 not in fine_gammas: fine_gammas.append(90)
            if min_g == 0 and 0 not in fine_gammas: fine_gammas.append(0)
            fine_gammas = sorted(list(set([round(g, 2) for g in fine_gammas if 0 <= g <= 90])))
            
            best_fine_torque = best_coarse_torque 
            
            for g in fine_gammas:
                if g in coarse_gammas: continue
                print(f" -> Fine: Gamma = {g}°")
                tq, res_arr = evaluate_gamma(g)
                if tq > best_fine_torque:
                    best_fine_torque = tq; best_res_array = res_arr
        else:
            best_idx = len(grid_map_data["Gamma_deg"]) - len(coarse_gammas) + coarse_gammas.index(best_coarse_gamma)
            best_res_array = [grid_map_data[k][best_idx] for k in master_mtpa_data.keys()]
                
        if best_res_array:
            for k, v in zip(master_mtpa_data.keys(), best_res_array):
                master_mtpa_data[k].append(v)

    # 3. Master Exports
    with open(os.path.join(mtpa_root, "Master_MTPA_Peaks.csv"), 'w', newline='') as f:
        w = csv.writer(f); w.writerow(master_mtpa_data.keys()); w.writerows(zip(*[master_mtpa_data[k] for k in master_mtpa_data.keys()]))
    with open(os.path.join(mtpa_root, "Master_Full_Grid_Map.csv"), 'w', newline='') as f:
        w = csv.writer(f); w.writerow(grid_map_data.keys()); w.writerows(zip(*[grid_map_data[k] for k in grid_map_data.keys()]))

    print("\n -> Compiling Massive 3D Database (.npz)...")
    unique_currents = np.unique([d["current"] for d in all_raw_data])
    unique_gammas = np.unique([d["gamma"] for d in all_raw_data])
    theta_axis = np.array(all_raw_data[0]["res"]["Theta_Elec_deg"])
    
    raw_torque_3d = np.full((len(unique_currents), len(unique_gammas), num_steps), np.nan)
    for d in all_raw_data:
        c_idx = np.where(unique_currents == d["current"])[0][0]
        g_idx = np.where(unique_gammas == d["gamma"])[0][0]
        raw_torque_3d[c_idx, g_idx, :] = d["res"]["Torque_Nm"]
        
    np.savez(
        os.path.join(mtpa_root, "Master_MTPA_Database.npz"),
        Axis_Current_A = unique_currents,
        Axis_Gamma_deg = unique_gammas,
        Axis_Theta_Elec_deg = theta_axis,
        **master_mtpa_data,
        Grid_I_d_A = grid_map_data["I_d_A"],
        Grid_I_q_A = grid_map_data["I_q_A"],
        Grid_Avg_Torque_Nm = grid_map_data["Avg_Torque_Nm"],
        Grid_V_d_Volts = grid_map_data["V_d_Volts"],
        Grid_V_q_Volts = grid_map_data["V_q_Volts"],
        Raw_Torque_3D = raw_torque_3d
    )

    print(f">>> MTPA Suite Complete! Master files and Database generated.")
    return master_mtpa_data