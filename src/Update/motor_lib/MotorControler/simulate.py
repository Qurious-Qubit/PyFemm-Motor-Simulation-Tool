import os
import math
import numpy as np

from .transient import execute_transient
from .mtpa import execute_mtpa

class MotorFEAEngine:
    def __init__(
        self, *, filename, base_folder, p, q, 
        initial_rotor_pos, IP, band_name, max_cores=None
    ):
        """
        Phase 1: Initialize the invariant 'Motor DNA'.
        band_name is now permanently stored in the DNA.
        """
        self.base_filename = filename
        self.base_folder = base_folder
        self.p = p
        self.q = q
        self.initial_rotor_pos = initial_rotor_pos
        self.IP = IP
        self.band_name = band_name
        
        self.max_cores = max_cores or max(1, os.cpu_count() - 1)
        self.BM = np.array([q, p]) / math.gcd(q, p)

    def _get_active_filename(self, name_extension):
        """Helper to stitch name extensions onto the base filename."""
        if not name_extension:
            return self.base_filename
        
        name_parts = os.path.splitext(self.base_filename)
        return f"{name_parts[0]}_{name_extension}{name_parts[1]}"

    # ==========================================
    # PHASE 2: EXECUTION APIS
    # ==========================================
    def run_transient(
        self, *, I_rms, gamma_elec_deg, theta_start_elec, theta_end_elec, num_steps, 
        Ns_rpm, delete_files, name_extension=None
    ):
        """Runs a single operating point transient analysis."""
        active_filename = self._get_active_filename(name_extension)
        
        return execute_transient(
            motor_state=self, 
            I_rms=I_rms, 
            gamma_elec_deg=gamma_elec_deg, 
            theta_start_elec=theta_start_elec, 
            theta_end_elec=theta_end_elec, 
            num_steps=num_steps, 
            Ns_rpm=Ns_rpm,
            active_filename=active_filename, 
            band_name=self.band_name,
            delete_files=delete_files
        )

    def run_mtpa(
        self, *, I_rms_list, theta_start_elec, theta_end_elec, num_steps, 
        Ns_rpm, delete_files,
        coarse_step_deg=15.0, run_fine_search=True, fine_window_deg=10.0, fine_step_deg=1.0,
        name_extension=None
    ):
        """Runs an intelligent MTPA sweep with protected boundaries."""
        active_filename = self._get_active_filename(name_extension)
        
        return execute_mtpa(
            motor_state=self, 
            I_rms_list=I_rms_list, 
            theta_start_elec=theta_start_elec, 
            theta_end_elec=theta_end_elec, 
            num_steps=num_steps, 
            Ns_rpm=Ns_rpm,
            coarse_step_deg=coarse_step_deg, 
            run_fine_search=run_fine_search, 
            fine_window_deg=fine_window_deg, 
            fine_step_deg=fine_step_deg,
            active_filename=active_filename, 
            band_name=self.band_name,
            delete_files=delete_files
        )