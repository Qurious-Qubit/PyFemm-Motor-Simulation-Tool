import os
import time
import math
import concurrent.futures
import numpy as np
import femm

def parallel_analyze_step(args):
    """Executes a single FEMM time step independently across CPU cores."""
    step_num, theta_mech_deg, i_a, i_b, i_c, p, BM, base_filename, step_filename, band_name, delete_files, max_cores = args
    
    stagger_delay = (step_num % max_cores) * 0.5
    time.sleep(stagger_delay)
    
    try:
        femm.openfemm(1)  
        femm.opendocument(base_filename)
        
        femm.mi_modifycircprop("fase1", 1, i_a)
        femm.mi_modifycircprop("fase2", 1, i_b)
        femm.mi_modifycircprop("fase3", 1, i_c)
        femm.mi_modifyboundprop(band_name, 10, theta_mech_deg)
        
        femm.mi_saveas(step_filename)
        femm.mi_analyze(1) 
        femm.mi_loadsolution()
        
        torque_total = femm.mo_gapintegral(band_name, 0)
        sym_mult = p / BM[1]
        
        props_a = femm.mo_getcircuitproperties("fase1")
        props_b = femm.mo_getcircuitproperties("fase2")
        props_c = femm.mo_getcircuitproperties("fase3")
        
        flux_a = props_a[2] * sym_mult
        flux_b = props_b[2] * sym_mult
        flux_c = props_c[2] * sym_mult
        
        femm.mo_close()
        femm.mi_close()
        femm.closefemm()
        
        if delete_files:
            try:
                if os.path.exists(step_filename): os.remove(step_filename)
                ans_filename = step_filename.replace('.fem', '.ans')
                if os.path.exists(ans_filename): os.remove(ans_filename)
            except: pass
                
        print(f"     [Core Task] Step {step_num:02d} complete. Torque: {torque_total:0.4f} N.m")
        return {
            "step_index": step_num, "Torque_Nm": torque_total,
            "I_A": props_a[0], "I_B": props_b[0], "I_C": props_c[0],
            "Flux_A_Wb": flux_a, "Flux_B_Wb": flux_b, "Flux_C_Wb": flux_c
        }
    except Exception as e:
        print(f"ERROR on step {step_num}: {e}")
        try: femm.closefemm() 
        except: pass
        return None

def post_process_dq0(*, results):
    """Converts ABC fluxes and currents to DQ0 frame."""
    theta_rad = np.radians(results["Theta_Elec_deg"])
    
    i_a = np.array(results["I_A"]); i_b = np.array(results["I_B"]); i_c = np.array(results["I_C"])
    f_a = np.array(results["Flux_A_Wb"]); f_b = np.array(results["Flux_B_Wb"]); f_c = np.array(results["Flux_C_Wb"])
    
    cos_t = np.cos(theta_rad); cos_t_120 = np.cos(theta_rad - 2*np.pi/3); cos_t_240 = np.cos(theta_rad + 2*np.pi/3)
    sin_t = np.sin(theta_rad); sin_t_120 = np.sin(theta_rad - 2*np.pi/3); sin_t_240 = np.sin(theta_rad + 2*np.pi/3)
    
    results["I_d_A"] = (2/3) * (i_a * cos_t + i_b * cos_t_120 + i_c * cos_t_240)
    results["I_q_A"] = -(2/3) * (i_a * sin_t + i_b * sin_t_120 + i_c * sin_t_240)
    
    results["Flux_d_Wb"] = (2/3) * (f_a * cos_t + f_b * cos_t_120 + f_c * cos_t_240)
    results["Flux_q_Wb"] = -(2/3) * (f_a * sin_t + f_b * sin_t_120 + f_c * sin_t_240)
        
    return results

def run_base_fea_sweep(
    *, motor_state, I_rms, gamma_elec_deg, theta_start_elec, theta_end_elec, num_steps, 
    active_filename, temp_folder, band_name, delete_files
):
    """Orchestrates the parallel sweep using parameters from the MotorFEAEngine."""
    results = {
        "Theta_Elec_deg": [], "Theta_Mech_deg": [], "Torque_Nm": [],
        "I_A": [], "I_B": [], "I_C": [], "Flux_A_Wb": [], "Flux_B_Wb": [], "Flux_C_Wb": []
    }
    
    I_peak = I_rms * math.sqrt(2)
    step_size_elec = (theta_end_elec - theta_start_elec) / (num_steps - 1) if num_steps > 1 else 0
    
    tasks = []
    for step in range(num_steps):
        theta_elec = theta_start_elec + (step * step_size_elec)
        theta_mech = theta_elec / (motor_state.p / 2) + motor_state.initial_rotor_pos
        
        i_a = I_peak * math.cos(math.radians(theta_elec + gamma_elec_deg + motor_state.IP))
        i_b = I_peak * math.cos(math.radians(theta_elec - 120 + gamma_elec_deg + motor_state.IP))
        i_c = I_peak * math.cos(math.radians(theta_elec + 120 + gamma_elec_deg + motor_state.IP))
        
        step_filename = os.path.join(temp_folder, f"step_{step:03d}.fem")
        tasks.append((step, theta_mech, i_a, i_b, i_c, motor_state.p, motor_state.BM, 
                      active_filename, step_filename, band_name, delete_files, motor_state.max_cores))
        
        results["Theta_Elec_deg"].append(theta_elec)
        results["Theta_Mech_deg"].append(theta_mech)

    with concurrent.futures.ProcessPoolExecutor(max_workers=motor_state.max_cores) as executor:
        unordered_results = list(executor.map(parallel_analyze_step, tasks))

    valid_results = [r for r in unordered_results if r is not None]
    for r in sorted(valid_results, key=lambda x: x["step_index"]):
        for k in ["Torque_Nm", "I_A", "I_B", "I_C", "Flux_A_Wb", "Flux_B_Wb", "Flux_C_Wb"]:
            results[k].append(r[k])

    return post_process_dq0(results=results)