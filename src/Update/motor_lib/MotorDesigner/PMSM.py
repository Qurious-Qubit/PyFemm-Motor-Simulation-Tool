import math

class PMSMDesigner:
    """
    Analytical design engine for Permanent Magnet Synchronous Machines (PMSM).
    Call PMSMDesigner.show_inputs() to see required and optional configurations.
    """

    # --- Self-Documenting Dictionaries ---
    _REQUIRED_INPUTS = {
        "Pm": "Mechanical power (kW)",
        "Ns": "Speed (RPM)",
        "B_avg": "Specific Magnetic Loading (T)",
        "ac": "Specific Electric Loading (A/mm)",
        "Vp": "Phase Voltage (V)",
        "q": "Number of slots",
        "p": "Number of poles",
        "Np": "Parallel paths",
        "bs0": "Slot opening width (mm)",
        "g": "Air-gap (mm)",
        "Br": "Remanent flux density of magnet (T)"
    }

    _OPTIONAL_INPUTS = {
        "mu_r": {"default": 1.05, "desc": "Relative permeability of magnet"},
        "Bm_target_ratio": {"default": 0.7, "desc": "Target operating point ratio (0.5 = Max Energy, 0.8 = Safe)"},
        "Jsw": {"default": 5.0, "desc": "Current density (A/mm^2)"},
        "kf": {"default": 0.35, "desc": "Copper Fill factor"},
        "ki": {"default": 0.95, "desc": "Iron Stacking factor"},
        "gamma_emf": {"default": 0.91, "desc": "Eph/Vph ratio"},
        "hs0": {"default": 0.3, "desc": "Slot opening depth (mm)"},
        "hs1": {"default": 0.5, "desc": "Slot opening wedge depth (mm)"},
        "Bst": {"default": 1.7, "desc": "Maximum tooth flux density (T)"},
        "Bsy": {"default": 1.4, "desc": "Maximum stator yoke flux density (T)"},
        "m": {"default": 3, "desc": "Number of phases"},
        "winding_layers": {"default": "DOUBLE", "desc": "Winding layers ('SINGLE' or 'DOUBLE')"},
        "n_eff": {"default": 0.92, "desc": "Efficiency"},
        "PF": {"default": 0.8, "desc": "Power Factor"},
        "kw": {"default": 0.975, "desc": "Winding Factor"},
        "ar_m": {"default": 0.888, "desc": "L/D ratio"}
    }

    @classmethod
    def show_inputs(cls):
        """Prints a beautifully formatted table of required and optional inputs."""
        print("="*60)
        print("=== PMSM DESIGNER: REQUIRED INPUTS ===")
        print("="*60)
        for key, desc in cls._REQUIRED_INPUTS.items():
            print(f"{key:<15} : {desc}")
        
        print("\n" + "="*60)
        print("=== PMSM DESIGNER: OPTIONAL INPUTS (Defaults shown) ===")
        print("="*60)
        for key, info in cls._OPTIONAL_INPUTS.items():
            print(f"{key:<15} : Default = {info['default']:<8} | {info['desc']}")
        print("="*60 + "\n")

    @classmethod
    def show_design_guidelines(cls):
        """Prints industry-standard design parameters for quick reference."""
        print("="*75)
        print("=== DESIGN GUIDELINES: SPECIFIC MAGNETIC LOADING (B_avg) ===")
        print("="*75)
        print("Induction Motor                   : 0.40 - 0.60 T")
        print("PMSM (Ferrite Magnets)            : 0.35 - 0.55 T")
        print("PMSM (NdFeB Surface Mount)        : 0.65 - 0.85 T")
        print("PMSM (NdFeB Interior / V-Type)    : 0.75 - 0.95 T")
        
        print("\n" + "="*75)
        print("=== DESIGN GUIDELINES: SPECIFIC ELECTRIC LOADING (ac) ===")
        print("="*75)
        print("Natural Convection (Air)          : 10 - 15 A/mm")
        print("Forced Fan Cooled (TEFC)          : 15 - 30 A/mm")
        print("Water/Liquid Jacket Cooling       : 30 - 60 A/mm")
        print("Direct Oil/Stator Flooding        : 60 - 90 A/mm")
        
        print("\n" + "="*75)
        print("=== DESIGN GUIDELINES: TOOTH FLUX DENSITY (Bst) ===")
        print("="*75)
        print("Typical Steel Knee Point          : ~ 1.50 - 1.60 T")
        print("Typical Hard Saturation Point     : ~ 2.00 - 2.10 T")
        print("Target Bst (Standard)             : 1.60 - 1.80 T")
        print("Target Bst (High-Perf/Traction)   : 1.80 - 2.00 T")
        print("-> Note: Teeth are intentionally pushed near or slightly above the")
        print("   knee point into mild saturation to maximize the copper slot area.")

        print("\n" + "="*75)
        print("=== DESIGN GUIDELINES: YOKE FLUX DENSITY (Bsy) ===")
        print("="*75)
        print("Target Bsy (Standard)             : 1.20 - 1.50 T")
        print("Target Bsy (High-Perf/Traction)   : 1.40 - 1.60 T")
        print("-> Note: The yoke is kept strictly below or at the knee point.")
        print("   Pushing the yoke into saturation drastically increases core losses.")

        print("\n" + "="*75)
        print("=== DESIGN GUIDELINES: CURRENT DENSITY (Jsw) ===")
        print("="*75)
        print("Natural Convection (Air)          : 3.0 - 5.0  A/mm²")
        print("Forced Fan Cooled (TEFC)          : 5.0 - 8.0  A/mm²")
        print("Water/Liquid Jacket Cooling       : 10.0 - 15.0 A/mm²")
        print("Direct Oil/Stator Flooding        : 15.0 - 30.0 A/mm²")
        print("="*75 + "\n")
        
    def __init__(
        self, 
        *, 
        # Required
        Pm, Ns, B_avg, ac, Vp, q, p, Np, bs0, g, Br,
        
        # Optional
        mu_r=1.05, Bm_target_ratio=0.7, Jsw=5.0, kf=0.35, ki=0.95, 
        gamma_emf=0.91, hs0=0.3, hs1=0.5, Bst=1.7, Bsy=1.4, m=3, 
        winding_layers="DOUBLE", n_eff=0.92, PF=0.8, kw=0.975, ar_m=0.888
    ):
        """Phase 1: Initialize the invariant 'Motor DNA'."""
        self.Pm = Pm
        self.Ns = Ns
        self.B_avg = B_avg
        self.ac = ac
        self.Vp = Vp
        self.q = q
        self.p = p
        self.Np = Np
        self.bs0 = bs0
        self.g = g
        self.Br = Br
        
        self.mu_r = mu_r
        self.Bm_target_ratio = Bm_target_ratio
        self.Jsw = Jsw
        self.kf = kf
        self.ki = ki
        self.gamma_emf = gamma_emf
        self.hs0 = hs0
        self.hs1 = hs1
        self.Bst = Bst
        self.Bsy = Bsy
        self.m = m
        self.winding_layers = winding_layers
        self.n_eff = n_eff
        self.PF = PF
        self.kw = kw
        self.ar_m = ar_m

        # Internal state to hold results
        self.sizing_results = None
        self.elec_results = None
        self.slot_results = None
        self.magnet_results = None

    def _machine_parameters(self):
        nrps = self.Ns / 60
        Pe = self.Pm / self.n_eff 
        P_g = Pe / self.PF 
        
        Trq = P_g / (2 * math.pi * nrps) * 1000

        G = 1.11 * (math.pi**2) * self.kw * self.B_avg * self.ac 
        D2L = (P_g / G / nrps) * 1e9  
        
        D = math.ceil((D2L / self.ar_m)**(1/3))
        L = math.ceil(D * self.ar_m)
        D2L_Actual = (D**2) * L

        Vr_m3 = (math.pi / 4) * D2L_Actual * 1e-9 
        TRV_kNm_m3 = (Trq / Vr_m3) / 1000 

        return {
            "D": D, "L": L, "G": G, "Pm": self.Pm, "Pe": Pe, 
            "Ns": self.Ns, "B_avg": self.B_avg, "kw": self.kw, 
            "D2L_Required": D2L, "D2L_Actual": D2L_Actual, 
            "PF": self.PF, "Torque": Trq, "TRV_kNm_m3": TRV_kNm_m3
        }

    def _electrical_parameters(self):
        D = self.sizing_results["D"]
        L = self.sizing_results["L"]
        Pe = self.sizing_results["Pe"]
        
        f = self.p * self.Ns / 120
        phi = self.B_avg * math.pi * D * L / 1e6 * 1000  
        phi_p = phi / self.p                             
        
        Ip = 1000 * Pe / self.PF / self.m / self.Vp 
        
        n_layers = 1 if self.winding_layers == "SINGLE" else 2 if self.winding_layers == "DOUBLE" else 0 
        Nc = (self.q / self.m / 2) if self.winding_layers == "SINGLE" else (self.q / self.m) if self.winding_layers == "DOUBLE" else 0 
        Nc_p = Nc / self.Np 
        Ic = Ip / self.Np 
        
        return {
            "phi_p": phi_p, "f": f, "Vp": self.Vp, "Ip": Ip, 
            "Ic": Ic, "Nc": Nc, "Nc_p": Nc_p, "n_layers": n_layers, 
            "q": self.q, "p": self.p, "m": self.m, "Np": self.Np
        }

    def _slot_design(self):
        D = self.sizing_results["D"]
        L = self.sizing_results["L"]
        
        phi_p = self.elec_results["phi_p"]
        f = self.elec_results["f"]
        Ic = self.elec_results["Ic"]
        Nc_p = self.elec_results["Nc_p"]
        n_layers = self.elec_results["n_layers"]

        phi_st_max = (phi_p / 2) * math.sin(math.pi * self.p / self.q) 
        wt = (phi_st_max / (self.Bst * L * self.ki)) * 1000 
        wsy = (phi_p / 2) / (L / 1000 * self.ki * self.Bsy)

        Ep = self.gamma_emf * self.Vp 
        Nt_ph = Ep / (4.44 * f * self.kw * phi_p) * 1000 
        Nt_c = Nt_ph / Nc_p 

        cAsc = Ic / self.Jsw              
        cAca = cAsc * Nt_c           
        gAca = cAca / self.kf             
        sA = n_layers * gAca         

        tan_pi_q = math.tan(math.pi / self.q)
        cos_pi_q = math.cos(math.pi / self.q)

        bs1 = 2 * ((tan_pi_q * (D / 2 + self.hs0 + self.hs1)) - ((wt / 2) / cos_pi_q))
        d = (-bs1 + math.sqrt((bs1**2) + (4 * tan_pi_q * sA))) / (2 * tan_pi_q)
        bs2 = bs1 + 2 * d * tan_pi_q

        return {
            "wt": wt, "wsy": wsy, "Ntc": Nt_c, "sA": sA, 
            "bs1": bs1, "bs2": bs2, "d": d, "bs0": self.bs0, "hs0": self.hs0
        }

    def _magnet_design(self):
        D = self.sizing_results["D"]
        L = self.sizing_results["L"]
        
        phi_p = self.elec_results["phi_p"]
        Ic = self.elec_results["Ic"]
        Nc = self.elec_results["Nc"]  
        Ntc = self.slot_results["Ntc"]

        u = self.bs0 / (2 * self.g)
        gamma = (4 / math.pi) * (u * math.atan(u) - math.log(math.sqrt(1 + u**2)))
        tau_s = (math.pi * D) / self.q
        kc = tau_s / (tau_s - gamma * self.g)
        g_effective = kc * self.g

        F_a = (self.m / 2) * (4 / math.pi) * (math.sqrt(2) * Ic * Ntc * Nc) / self.p

        Bm_target = self.Br * self.Bm_target_ratio
        
        phi_webers = phi_p / 1000 
        wm_total_required = (phi_webers / Bm_target) * 1e6 / (L * self.ki)

        mu_0 = 4 * math.pi * 1e-7
        g_meters = g_effective / 1000
        B_peak_gap = self.B_avg * (math.pi / 2)
        
        hm_meters = (self.mu_r * ((B_peak_gap * g_meters / mu_0) + F_a)) / ((self.Br - Bm_target) / mu_0)
        hm_opt = hm_meters * 1000

        tau_p = (math.pi * D) / self.p
        max_physical_v_width = 1.8 * tau_p
        
        design_feasible = True
        warning_msg = "OK"

        if wm_total_required > max_physical_v_width:
            design_feasible = False
            warning_msg = (
                f"WARNING: Required magnet width ({wm_total_required:.1f} mm) exceeds "
                f"estimated physical space ({max_physical_v_width:.1f} mm)."
            )

        return {
            "kc": kc, "Fa": F_a, "Bm_Target": Bm_target, 
            "wm_opt": wm_total_required, "hm_opt": hm_opt, "Lm": L, 
            "Design_Feasible_Flag": design_feasible, "Warning_Message": warning_msg
        }

    def run_design(self):
        """Phase 2: Executes the full pipeline sequentially based on the DNA."""
        print("=== EXECUTING MOTOR ANALYTICAL DESIGN ===\n")
        
        self.sizing_results = self._machine_parameters()
        self.elec_results = self._electrical_parameters()
        self.slot_results = self._slot_design()
        self.magnet_results = self._magnet_design()

        return {
            "sizing_results": self.sizing_results,
            "electrical_results": self.elec_results,
            "slot_results": self.slot_results,
            "magnet_results": self.magnet_results
        }


# ==========================================
# EXAMPLE EXECUTION (Can be run anywhere)
# ==========================================
if __name__ == "__main__":
    
    # 1. Self-Documenting Check: Find out what we need before coding!
    PMSMDesigner.show_inputs()

    # 2. Initialize the Class
    my_designer = PMSMDesigner(
        Pm=10.0,
        Ns=10000,
        B_avg=(2.0 / math.pi) * 0.9,
        ac=30,
        Vp=230.0,
        q=48,
        p=8,
        Np=2,
        bs0=0.2,
        g=0.8,
        Br=1.2,
        Bm_target_ratio=0.5  # Overriding an optional parameter
    )

    # 3. Run the design pipeline
    results = my_designer.run_design()
    
    # 4. Access output
    print(f"Calculated Air-Gap Diameter (D) : {results['sizing_results']['D']} mm")
    print(f"Calculated Stack Length (L)     : {results['sizing_results']['L']} mm")
    print(f"Feasibility Status              : {results['magnet_results']['Warning_Message']}")