import os
import sys

# Point Python directly to your Dev folder so it can find 'motor_lib'
sys.path.append(r"C:\D_Drive\Technical\MyPapers\CoggingVsLeakage\ModuleLibr\Dev")

from motor_lib.MotorControler.simulate import MotorFEAEngine

def run_strict_external_test():
    current_dir = os.getcwd()
    base_filename = "IPM_Motor_Base.fem" 

    print("==================================================")
    print(f"TESTING STRICT MOTOR FEA ENGINE")
    print(f"Output Directory: {current_dir}")
    print("==================================================")

    # 1. Initialize the Motor Object (The DNA)
    my_motor = MotorFEAEngine(
        filename=base_filename,
        base_folder=current_dir,
        p=8,
        q=48,
        initial_rotor_pos=3.75,
        IP=90.0,
        band_name="AGap"
    )
    
    # 2. Run a Quick Transient Simulation
    print("\n>>> Executing Test: TRANSIENT")
    my_motor.run_transient(
        I_rms=10.0,
        gamma_elec_deg=0.0,
        theta_start_elec=0,
        theta_end_elec=60,
        num_steps=11,        
        Ns_rpm=10000,        
        delete_files=True
    )
    
    # 3. Run a Quick MTPA Search
    print("\n>>> Executing Test: MTPA")
    my_motor.run_mtpa(
        I_rms_list=[10.0, 20.0],
        theta_start_elec=0,
        theta_end_elec=60,
        num_steps=11,        
        Ns_rpm=10000,        
        delete_files=True,
        coarse_step_deg=15.0,
        run_fine_search=True,
        fine_window_deg=10.0,
        fine_step_deg=2.0    
    )
    
    print("\nAll external tests completed successfully! Check the generated folders.")

if __name__ == '__main__':
    run_strict_external_test()