import numpy as np
import math
import femm

from .femmStarter import startBuild
from .buildStator import draw_stator
from .buildRotor import draw_rotor

class VTypeBuilder:
    """
    Geometry builder for V-Type Interior Permanent Magnet (IPM) Motors.
    Call VTypeBuilder.show_inputs() to see all required and optional parameters.
    """

    # --- Self-Documenting Dictionaries ---
    _REQUIRED_INPUTS = {
        "filename": "Active FEMM file name (e.g., 'motor.fem')",
        "L": "Stack length (mm)",
        "m": "Number of phases",
        "D": "Stator bore diameter (mm)",
        "Ds": "Stator outer diameter (mm)",
        "p": "Number of poles",
        "q": "Number of slots",
        "g": "Air-gap size (mm)",
        "pole_arc_elec_deg": "Pole arc in electrical degrees",
        "mag_thickness": "Magnet thickness (mm)",
        "mag_width": "Magnet width per pole (mm)",
        "mag_pkt_width": "Magnet pocket width (mm)",
        "tr_t": "Rotor bridge thickness (mm)",
        "hs0": "Slot opening depth (mm)",
        "hs1": "Slot wedge depth (mm)",
        "d": "Slot depth (mm)",
        "bs0": "Slot opening width (mm)",
        "bs1": "Slot top width (mm)",
        "bs2": "Slot bottom width (mm)",
        "wsy": "Stator yoke thickness (mm)",
        "winding_layout": "2D Numpy array of the winding layout",
        "Nt_c": "Turns per coil",
        "stator_iron_mat": "Name of stator iron material",
        "coil_mat": "Name of coil material",
        "rotor_iron_mat": "Name of rotor iron material",
        "magnet_mat": "Name of magnet material"
    }

    _OPTIONAL_INPUTS = {
        "standard_materials": {"default": "None", "desc": "List of standard FEMM materials to import"},
        "custom_bh_materials": {"default": "None", "desc": "Dict of custom BH cores: {'Name': 'path.csv'}"},
        "custom_magnets": {"default": "None", "desc": "Dict of custom magnets: {'Name': {'ur': 1.05, 'H_c': 490000}}"},
        "custom_conductors": {"default": "None", "desc": "Dict of custom conductors: {'Name': Cduct}"},
        "agap_maxsegdeg": {"default": 1.0, "desc": "Max degrees per arc segment for airgap sliding band"},
        "mesh_agap": {"default": 0.15, "desc": "Mesh size for airgap"},
        "mesh_stator": {"default": "None", "desc": "Mesh size for stator core"},
        "mesh_coil": {"default": "None", "desc": "Mesh size for coils"},
        "mesh_rotor_core": {"default": 1.0, "desc": "Mesh size for inner rotor core"},
        "mesh_rotor_embedded_pole": {"default": 0.75, "desc": "Mesh size for rotor pole pieces"},
        "mesh_rib": {"default": 0.15, "desc": "Mesh size for rotor iron bridges"},
        "mesh_magnet": {"default": 0.5, "desc": "Mesh size for magnets"},
        "mesh_shaft": {"default": 5.0, "desc": "Mesh size for shaft air"}
    }

    @classmethod
    def show_inputs(cls):
        """Prints a beautifully formatted table of required and optional inputs."""
        print("="*70)
        print("=== V-TYPE GEOMETRY BUILDER: REQUIRED INPUTS ===")
        print("="*70)
        for key, desc in cls._REQUIRED_INPUTS.items():
            print(f"{key:<20} : {desc}")
        
        print("\n" + "="*70)
        print("=== V-TYPE GEOMETRY BUILDER: OPTIONAL INPUTS (Defaults shown) ===")
        print("="*70)
        for key, info in cls._OPTIONAL_INPUTS.items():
            print(f"{key:<25} : Default = {str(info['default']):<6} | {info['desc']}")
        print("="*70 + "\n")

    def __init__(
        self,
        *,
        # --- Required Base Geometry & Parameters ---
        filename, L, m, D, Ds, p, q, g, 
        pole_arc_elec_deg, mag_thickness, mag_width, mag_pkt_width, tr_t,
        hs0, hs1, d, bs0, bs1, bs2, wsy,
        winding_layout, Nt_c,
        stator_iron_mat, coil_mat, rotor_iron_mat, magnet_mat,

        # --- Optional Material Dictionaries ---
        standard_materials=None,
        custom_bh_materials=None,
        custom_magnets=None,
        custom_conductors=None,

        # --- Optional Meshing Parameters ---
        agap_maxsegdeg=1.0,
        mesh_agap=0.15,
        mesh_stator=None,
        mesh_coil=None,
        mesh_rotor_core=1.0,
        mesh_rotor_embedded_pole=0.75,
        mesh_rib=0.15,
        mesh_magnet=0.5,
        mesh_shaft=5.0
    ):
        """Phase 1: Initialize the invariant 'Motor DNA'."""
        self.filename = filename
        self.L = L
        self.m = m
        self.D = D
        self.Ds = Ds
        self.p = p
        self.q = q
        self.g = g
        self.pole_arc_elec_deg = pole_arc_elec_deg
        self.mag_thickness = mag_thickness
        self.mag_width = mag_width
        self.mag_pkt_width = mag_pkt_width
        self.tr_t = tr_t
        self.hs0 = hs0
        self.hs1 = hs1
        self.d = d
        self.bs0 = bs0
        self.bs1 = bs1
        self.bs2 = bs2
        self.wsy = wsy
        self.winding_layout = winding_layout
        self.Nt_c = Nt_c
        
        self.stator_iron_mat = stator_iron_mat
        self.coil_mat = coil_mat
        self.rotor_iron_mat = rotor_iron_mat
        self.magnet_mat = magnet_mat

        self.standard_materials = standard_materials
        self.custom_bh_materials = custom_bh_materials
        self.custom_magnets = custom_magnets
        self.custom_conductors = custom_conductors

        self.agap_maxsegdeg = agap_maxsegdeg
        self.mesh_agap = mesh_agap
        self.mesh_stator = mesh_stator
        self.mesh_coil = mesh_coil
        self.mesh_rotor_core = mesh_rotor_core
        self.mesh_rotor_embedded_pole = mesh_rotor_embedded_pole
        self.mesh_rib = mesh_rib
        self.mesh_magnet = mesh_magnet
        self.mesh_shaft = mesh_shaft

        # Derived Parameters calculated once and stored
        self.Dr = self.D - (2 * self.g)
        self.BM = np.array([self.q, self.p]) / math.gcd(self.q, self.p)

    def start_workspace(self):
        """Initializes the FEMM document and imports materials."""
        print(f"Initializing FEMM Workspace: {self.filename}...")
        startBuild(
            filename=self.filename, 
            L=self.L, 
            m=self.m, 
            standard_materials=self.standard_materials,
            custom_bh_materials=self.custom_bh_materials,
            custom_magnets=self.custom_magnets,
            custom_conductors=self.custom_conductors
        )

    def draw_stator_geometry(self):
        """Draws the stator based on the stored DNA."""
        print("Drawing Stator...")
        draw_stator(
            D=self.D, Ds=self.Ds, p=self.p, q=self.q, g=self.g, 
            hs0=self.hs0, hs1=self.hs1, d=self.d, bs0=self.bs0, bs1=self.bs1, bs2=self.bs2, wsy=self.wsy,
            BM=self.BM, winding=self.winding_layout, Nt_c=self.Nt_c, filename=self.filename,
            stator_iron_mat=self.stator_iron_mat, coil_mat=self.coil_mat,
            mesh_agap=self.mesh_agap, mesh_stator=self.mesh_stator, mesh_coil=self.mesh_coil
        )

    def draw_rotor_geometry(self):
        """Draws the V-type rotor based on the stored DNA."""
        print("Drawing Rotor...")
        draw_rotor(
            Dr=self.Dr, g=self.g, p=self.p, q=self.q, BM=self.BM, 
            pole_arc_elec_deg=self.pole_arc_elec_deg, mag_thickness=self.mag_thickness, 
            mag_width=self.mag_width, mag_pkt_width=self.mag_pkt_width, tr_t=self.tr_t,
            filename=self.filename,
            rotor_iron_mat=self.rotor_iron_mat, magnet_mat=self.magnet_mat,
            mesh_agap=self.mesh_agap, mesh_rotor_core=self.mesh_rotor_core, 
            mesh_rotor_embedded_pole=self.mesh_rotor_embedded_pole, 
            mesh_rib=self.mesh_rib, mesh_magnet=self.mesh_magnet, mesh_shaft=self.mesh_shaft
        )

    def apply_airgap_boundaries(self):
        """Binds the anti-periodic or periodic boundaries along the airgap."""
        print("Binding Airgap Elements...")
        femm.openfemm()
        femm.opendocument(self.filename)

        bore_radius = self.D / 2
        rotor_radius = bore_radius - self.g
        
        R_stator_agap = bore_radius - (self.g / 3)
        R_rotor_agap = rotor_radius + (self.g / 3)
        
        inner_angle_deg = int(self.BM[1]) * (360 / self.p)
        outer_angle_deg = int(self.BM[0]) * (360 / self.q)

        if int(self.BM[1]) < self.p:
            is_anti_periodic = (int(self.BM[1]) % 2 != 0)
            agap_type = 7 if is_anti_periodic else 6 
            
            try:
                femm.mi_addboundprop("AGap", 0, 0, 0, 0, 0, 0, 0, 0, agap_type, 0, 0) 
            except:
                pass

            mid_s_angle = math.radians(outer_angle_deg / 2)
            femm.mi_selectarcsegment(R_stator_agap * math.cos(mid_s_angle), R_stator_agap * math.sin(mid_s_angle))
            femm.mi_setarcsegmentprop(self.agap_maxsegdeg, "AGap", 0, 0)
            femm.mi_clearselected()

            mid_r_angle = math.radians(inner_angle_deg / 2)
            femm.mi_selectarcsegment(R_rotor_agap * math.cos(mid_r_angle), R_rotor_agap * math.sin(mid_r_angle))
            femm.mi_setarcsegmentprop(self.agap_maxsegdeg, "AGap", 0, 0)
            femm.mi_clearselected()
        else:
            print("WARNING (AGap): Full machine modeled. AGap boundaries skipped.")

        femm.mi_saveas(self.filename)
        femm.mi_close()
        femm.closefemm()

    def build_motor(self):
        """Executes the entire drawing pipeline in the correct order."""
        self.start_workspace()
        self.draw_stator_geometry()
        self.draw_rotor_geometry()
        self.apply_airgap_boundaries()
        print("\n>>> Simulation Geometry Build Complete! <<<")
        return True


# ==========================================
# USAGE EXAMPLE (Can be run anywhere)
# ==========================================
if __name__ == "__main__":
    
    # 1. Ask the class what it needs (Self-Documenting)
    VTypeBuilder.show_inputs()

    # 2. Setup your winding matrix
    winding_matrix = np.array([
        [1, 1, -3, -3, 2, 2],
        [1, -3, -3, 2, 2, -1]
    ])
    
    # 3. Initialize the Builder Object
    builder = VTypeBuilder(
        filename="IPM_Motor_Base.fem",
        L=150.0,
        m=3,
        D=100.0,
        Ds=150.0,
        p=8,
        q=48,
        g=0.8,
        pole_arc_elec_deg=120.0,
        mag_thickness=3.5,
        mag_width=12.0,
        mag_pkt_width=20.0,
        tr_t=1.8,
        hs0=0.5,
        hs1=1.0,
        d=15.0,
        bs0=2.0,
        bs1=4.5,
        bs2=6.5,
        wsy=10.0,
        winding_layout=winding_matrix,
        Nt_c=12,
        stator_iron_mat="US Steel Type 2-S 0.024 inch thickness",
        coil_mat="24 SWG",
        rotor_iron_mat="US Steel Type 2-S 0.024 inch thickness",
        magnet_mat="N35",
        
        # Adding optional materials
        standard_materials=["Air", "24 SWG", "N35", "US Steel Type 2-S 0.024 inch thickness"]
    )

    # 4. Execute the build (You can call these step-by-step for debugging, or just use build_motor())
    builder.build_motor()