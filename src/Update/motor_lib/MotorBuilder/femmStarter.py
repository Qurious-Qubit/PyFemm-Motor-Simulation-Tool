import femm
import os
import csv

def addCustomBH(
    *,
    mat_name, 
    file_path, 
    lam_fill=1
):
    
    mu_x=1         # Relative permeability in x/r (Defaults to 1 for BH curve)
    mu_y=1         # Relative permeability in y/z
    H_c=0          # Coercivity (A/m)
    J=0            # Source current density (MA/m^2)
    Cduct=0        # Electrical conductivity (MS/m)
    Lam_d=0        # Lamination thickness (mm)
    Phi_hmax=0     # Hysteresis lag angle
    LamType=0      # Lamination type (0=None, 1=parallel to x, 2=parallel to y)
    Phi_hx=0,       # Hysteresis lag in x
    Phi_hy=0        # Hysteresis lag in y

    if not os.path.exists(file_path):
        print(f"Error: Custom material file '{file_path}' not found.")
        return False
        
    B_array = []
    H_array = []
    
    # Parse the file
    with open(file_path, 'r') as f:
        if file_path.endswith('.csv'):
            reader = csv.reader(f)
            for row in reader:
                try:
                    B_array.append(float(row[0]))
                    H_array.append(float(row[1]))
                except ValueError:
                    continue
        else:
            for line in f:
                parts = line.strip().split()
                if len(parts) >= 2:
                    try:
                        B_array.append(float(parts[0].replace(',', '')))
                        H_array.append(float(parts[1].replace(',', '')))
                    except ValueError:
                        continue

    if not B_array or not H_array:
        print(f"Error: Could not extract valid B-H data from '{file_path}'.")
        return False

    # 1. Initialize the custom material base with all passed parameters
    femm.mi_addmaterial(
        mat_name, mu_x, mu_y, H_c, J, Cduct, 
        Lam_d, Phi_hmax, lam_fill, LamType, Phi_hx, Phi_hy
    )
    
    # 2. Add B-H data points
    for b, h in zip(B_array, H_array):
        femm.mi_addbhpoint(mat_name, b, h)
        
    print(f"Successfully added Custom Core Material '{mat_name}' with {len(B_array)} B-H points.")
    return True

def addCustomMag(mat_name, ur, H_c):
    J=0            # Source current density (MA/m^2)
    Cduct=0        # Electrical conductivity (MS/m)
    Lam_d=0        # Lamination thickness (mm)
    Phi_hmax=0     # Hysteresis lag angle
    LamType=0      # Lamination type (0=None, 1=parallel to x, 2=parallel to y)
    Phi_hx=0       # Hysteresis lag in x
    Phi_hy=0        # Hysteresis lag in y
    lam_fill = 0    # Not relevant for magnets, but required by the function signature

    femm.mi_addmaterial(
        mat_name, ur, ur, H_c, J, Cduct, 
        Lam_d, Phi_hmax, lam_fill, LamType, Phi_hx, Phi_hy
    )

def addCustomCu(mat_name, Cduct):

    mu_x=1         # Relative permeability in x/r (Defaults to 1 for BH curve)
    mu_y=1         # Relative permeability in y/z
    H_c=0          # Coercivity (A/m)
    J=0            # Source current density (MA/m^2)
    Lam_d=0        # Lamination thickness (mm)
    Phi_hmax=0     # Hysteresis lag angle
    LamType=0      # Lamination type (0=None, 1=parallel to x, 2=parallel to y)
    Phi_hx=0       # Hysteresis lag in x
    Phi_hy=0        # Hysteresis lag in y
    lam_fill = 0    # Not relevant for coil conductors, but required by the function signature

    femm.mi_addmaterial(
        mat_name, mu_x, mu_y, H_c, J, Cduct, 
        Lam_d, Phi_hmax, lam_fill, LamType, Phi_hx, Phi_hy
    )
    return True

def startBuild(
    *, 
    filename, 
    L, 
    m, 
    standard_materials, 
    custom_bh_materials,  # dict: {"Core_Mat_Name": "path/to/bh.csv"}
    custom_magnets,       # dict: {"Mag_Name": {"ur": 1.05, "H_c": 490000}}
    custom_conductors        # dict: {"Cu_Name": Cduct_value}
):
    """
    Initializes a new FEMM magnetics document, defines the problem, 
    imports standard/custom materials, sets up phase circuits, and saves the file.
    """
    try:
        # Start FEMM and create a new Magnetics document
        femm.openfemm()
        femm.newdocument(0)
        
        # Problem definition: Frequency (0 for DC), Units, Type, Precision, Depth (L), Min Angle
        femm.mi_probdef(0, "millimeters", "planar", 1e-8, L, 30)

        # Disable smart mesh
        femm.mi_smartmesh(0)

        # --- 1. Import Standard Library Materials ---
        if standard_materials:
            for mat in standard_materials:
                femm.mi_getmaterial(mat)
                print(f"Imported standard material: {mat}")

        # --- 2. Import Custom B-H Non-Linear Materials ---
        if custom_bh_materials:
            for mat_name, file_path in custom_bh_materials.items():
                addCustomBH(mat_name=mat_name, file_path=file_path)

        # --- 3. Import Custom Magnets ---
        if custom_magnets:
            for mat_name, props in custom_magnets.items():
                # Expects a dict with 'ur' and 'H_c' keys for safety
                addCustomMag(
                    mat_name=mat_name, 
                    ur=props["ur"], 
                    H_c=props["H_c"]
                )

        # --- 4. Import Custom Copper / Conductors ---
        if custom_conductors:
            for mat_name, cduct in custom_conductors.items():
                addCustomCu(
                    mat_name=mat_name, 
                    Cduct=cduct
                )

        # --- 5. Create Phase Circuits ---
        for i in range(m):
            circuit_name = f"fase{i+1}"
            femm.mi_addcircprop(circuit_name, 0, 1) 
        print(f"Created {m} phase circuits.")

        # Save the file
        femm.mi_saveas(filename)
        print(f"FEMM initialization complete. Saved as: {filename}")
        
        return 1

    except Exception as e:
        print(f"An error occurred while initializing FEMM: {e}")
        return 0

