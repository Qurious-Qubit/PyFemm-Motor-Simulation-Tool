import femm
import math



def get_circle_line_intersection(R, m, c):
    A = 1 + m**2; B = 2 * m * c; C = c**2 - R**2
    det = B**2 - 4*A*C
    if det < 0: return 0, 0
    x_val = max((-B + math.sqrt(det)) / (2*A), (-B - math.sqrt(det)) / (2*A))
    return x_val, m * x_val + c

def draw_rotor(
    *,
    Dr, g, p, q, BM, pole_arc_elec_deg, mag_thickness, mag_width, mag_pkt_width, tr_t, filename, 
    rotor_iron_mat, magnet_mat, mesh_agap, mesh_rotor_core, mesh_rotor_embedded_pole, 
    mesh_rib, mesh_magnet, mesh_shaft
):
    femm.openfemm()
    femm.opendocument(filename)
    
    rotor_radius = Dr / 2
    shaft_radius = rotor_radius * 0.2 
    r_air = rotor_radius + (g / 3)
    
    theta_pole = 360 / p
    ang_pole_rad = math.radians(theta_pole / 2)
    pole_arc_mech_deg = pole_arc_elec_deg / (p / 2)
    theta_arc_rad = math.radians(pole_arc_mech_deg / 2)
    
    ref_design_radius = rotor_radius - tr_t
    
    x1_ref = rotor_radius * math.cos(theta_arc_rad)
    y1_ref = rotor_radius * math.sin(theta_arc_rad)
    
    min_len = ref_design_radius * math.sin(theta_arc_rad)
    pocket_length = mag_pkt_width
    if pocket_length <= min_len: pocket_length = y1_ref * 1.01 
        
    beta_rad = math.asin(y1_ref / pocket_length)
    m1 = math.tan(beta_rad) 
    
    x2 = x1_ref - (y1_ref / m1)
    shift_distance = mag_thickness / math.sin(beta_rad)
    x3 = x2 - shift_distance
    
    x_mid_in = (x2 + x3) / 2; y_mid_in = 0
    c_mid = -m1 * x_mid_in
    x_mid_out, y_mid_out = get_circle_line_intersection(rotor_radius, m1, c_mid)
    
    xc = (x_mid_in + x_mid_out) / 2; yc = (y_mid_out) / 2
    dx = x_mid_out - x_mid_in; dy = y_mid_out - y_mid_in
    actual_pocket_len = math.hypot(dx, dy)
    
    ux = dx / actual_pocket_len; uy = dy / actual_pocket_len
    nx = -uy; ny = ux   
    if (nx*xc + ny*yc) < 0: nx = uy; ny = -ux
        
    hl = mag_width / 2; ht = mag_thickness / 2
    mx1 = xc + hl*ux + ht*nx; my1 = yc + hl*uy + ht*ny 
    mx2 = xc - hl*ux + ht*nx; my2 = yc - hl*uy + ht*ny 
    mx3 = xc - hl*ux - ht*nx; my3 = yc - hl*uy - ht*ny 
    mx4 = xc + hl*ux - ht*nx; my4 = yc + hl*uy - ht*ny 
    
    design_radius = rotor_radius - tr_t
    c_outer = -m1 * x2
    x1, y1 = get_circle_line_intersection(design_radius, m1, c_outer)
    c_inner = -m1 * x3
    x4, y4 = get_circle_line_intersection(design_radius, m1, c_inner)
    
    cos_p = math.cos(ang_pole_rad); sin_p = math.sin(ang_pole_rad)
    x_sh_p = shaft_radius * cos_p;  y_sh_p = shaft_radius * sin_p
    x_des_p = design_radius * cos_p; y_des_p = design_radius * sin_p
    x_out_p = rotor_radius * cos_p;  y_out_p = rotor_radius * sin_p

    femm.mi_addnode(0, 0)
    femm.mi_addnode(x_sh_p, y_sh_p); femm.mi_addnode(x_des_p, y_des_p); femm.mi_addnode(x_out_p, y_out_p)
    femm.mi_addnode(x_sh_p, -y_sh_p); femm.mi_addnode(x_des_p, -y_des_p); femm.mi_addnode(x_out_p, -y_out_p)
    femm.mi_addsegment(0, 0, x_sh_p, y_sh_p); femm.mi_addsegment(x_sh_p, y_sh_p, x_des_p, y_des_p)
    femm.mi_addsegment(x_des_p, y_des_p, x_out_p, y_out_p)
    femm.mi_addsegment(0, 0, x_sh_p, -y_sh_p); femm.mi_addsegment(x_sh_p, -y_sh_p, x_des_p, -y_des_p)
    femm.mi_addsegment(x_des_p, -y_des_p, x_out_p, -y_out_p)

    femm.mi_addnode(x1, y1); femm.mi_addnode(x4, y4); femm.mi_addnode(x1, -y1); femm.mi_addnode(x4, -y4)
    femm.mi_addnode(x2, 0); femm.mi_addnode(x3, 0)

    for flip_y in [1, -1]:
        mx1_f = mx1; my1_f = my1 * flip_y; mx2_f = mx2; my2_f = my2 * flip_y
        mx3_f = mx3; my3_f = my3 * flip_y; mx4_f = mx4; my4_f = my4 * flip_y
        
        femm.mi_addnode(mx1_f, my1_f); femm.mi_addnode(mx2_f, my2_f)
        femm.mi_addnode(mx3_f, my3_f); femm.mi_addnode(mx4_f, my4_f)
        
        # Magnet Rectangle
        femm.mi_addsegment(mx1_f, my1_f, mx2_f, my2_f); femm.mi_addsegment(mx2_f, my2_f, mx3_f, my3_f)
        femm.mi_addsegment(mx3_f, my3_f, mx4_f, my4_f); femm.mi_addsegment(mx4_f, my4_f, mx1_f, my1_f)

        # --- THE FIX: Discrete segments to connect the pocket walls to the magnet ---
        y1_f = y1 * flip_y
        y4_f = y4 * flip_y
        
        # Outer pocket segments (Rib to Magnet Top)
        femm.mi_addsegment(x1, y1_f, mx1_f, my1_f)
        femm.mi_addsegment(x4, y4_f, mx4_f, my4_f)
        
        # Inner pocket segments (Magnet Bottom to Center Bridge)
        femm.mi_addsegment(mx2_f, my2_f, x2, 0)
        femm.mi_addsegment(mx3_f, my3_f, x3, 0)

    # (Removed the overlapping continuous full-length pocket segments from here)

    femm.mi_addarc(x_sh_p, -y_sh_p, x_sh_p, y_sh_p, theta_pole, 1)
    femm.mi_addarc(x_out_p, -y_out_p, x_out_p, y_out_p, theta_pole, 1)
    femm.mi_addarc(x_des_p, -y_des_p, x4, -y4, math.degrees(ang_pole_rad - math.atan2(y4, x4)), 1)
    femm.mi_addarc(x4, -y4, x1, -y1, math.degrees(math.atan2(y4, x4) - math.atan2(y1, x1)), 1)
    femm.mi_addarc(x1, -y1, x1, y1, math.degrees(2 * math.atan2(y1, x1)), 1)
    femm.mi_addarc(x1, y1, x4, y4, math.degrees(math.atan2(y4, x4) - math.atan2(y1, x1)), 1)
    femm.mi_addarc(x4, y4, x_des_p, y_des_p, math.degrees(ang_pole_rad - math.atan2(y4, x4)), 1)

    safe_selection_radius = rotor_radius + (g / 10) 
    
    femm.mi_selectcircle(0, 0, safe_selection_radius, 4) 
    femm.mi_moverotate(0, 0, theta_pole / 2)
    femm.mi_clearselected()
    
    if int(BM[1] - 1) > 0:
        femm.mi_selectcircle(0, 0, safe_selection_radius, 4) 
        femm.mi_copyrotate2(0, 0, theta_pole, int(BM[1] - 1), 4) 
        femm.mi_clearselected()

    # --- DRAW AIRGAP SLIDING BAND (Single Arc) ---
    total_angle_deg = int(BM[1]) * theta_pole
    
    if int(BM[1]) < p:
        xb = r_air * math.cos(math.radians(total_angle_deg))
        yb = r_air * math.sin(math.radians(total_angle_deg))
        femm.mi_addnode(r_air, 0); femm.mi_addnode(xb, yb)
        femm.mi_addarc(r_air, 0, xb, yb, total_angle_deg, 1)

        femm.mi_addsegment(rotor_radius, 0, r_air, 0)
        xe_out = rotor_radius * math.cos(math.radians(total_angle_deg))
        ye_out = rotor_radius * math.sin(math.radians(total_angle_deg))
        femm.mi_addsegment(xe_out, ye_out, xb, yb)
    else:
        femm.mi_addnode(r_air, 0); femm.mi_addnode(-r_air, 0)
        femm.mi_addarc(r_air, 0, -r_air, 0, 180, 1)
        femm.mi_addarc(-r_air, 0, r_air, 0, 180, 1)

    # --- BOUNDARIES (PERIODIC ONLY) ---
    is_anti_periodic = (int(BM[1]) % 2 != 0)
    b_type = 5 if is_anti_periodic else 4 
    
    try:
        if int(BM[1]) < p:
            for i in range(1, 5): femm.mi_addboundprop(f"r{i}", 0, 0, 0, 0, 0, 0, 0, 0, b_type)
    except:
        pass

    if int(BM[1]) < p:
        r_mids = [shaft_radius / 2, (shaft_radius + design_radius) / 2, (design_radius + rotor_radius) / 2, (rotor_radius + r_air) / 2]
        total_r_rad = math.radians(total_angle_deg)
        
        for i, r_mid in enumerate(r_mids):
            b_name = f"r{i+1}"
            femm.mi_selectsegment(r_mid, 0)
            femm.mi_setsegmentprop(b_name, 0, 1, 0, 0)
            femm.mi_clearselected()
            
            femm.mi_selectsegment(r_mid * math.cos(total_r_rad), r_mid * math.sin(total_r_rad))
            femm.mi_setsegmentprop(b_name, 0, 1, 0, 0)
            femm.mi_clearselected()
    else:
        print("WARNING (Rotor): Full machine modeled. Cyclic boundaries skipped.")

    # --- MATERIAL ASSIGNMENT ---
    mag_angle_base = math.degrees(math.atan2(ny, nx))
    lbl_shaft_air = (shaft_radius / 2, 0); lbl_inner_core = ((shaft_radius + x3) / 2, 0)
    lbl_central_air = ((x3 + x2) / 2, 0); lbl_pole_core = ((x2 + design_radius) / 2, 0)
    lbl_rib_core = ((design_radius + rotor_radius) / 2, 0)
    ax_out = (x_mid_out + (xc + hl*ux)) / 2; ay_out = (y_mid_out + (yc + hl*uy)) / 2
    
    lbl_airgap_x = ((rotor_radius + r_air) / 2) * math.cos(math.radians(theta_pole / 2))
    lbl_airgap_y = ((rotor_radius + r_air) / 2) * math.sin(math.radians(theta_pole / 2))
    femm.mi_addblocklabel(lbl_airgap_x, lbl_airgap_y)
    femm.mi_selectlabel(lbl_airgap_x, lbl_airgap_y)
    femm.mi_setblockprop("Air", 0, mesh_agap, "<None>", 0, 0, 0)
    femm.mi_clearselected()
    
    for i in range(int(BM[1])):
        rot_angle = math.radians((theta_pole / 2) + (i * theta_pole))
        cos_a = math.cos(rot_angle); sin_a = math.sin(rot_angle)
        
        def place_label(x, y, mat, mesh_size, prop_angle=0):
            rx = x * cos_a - y * sin_a; ry = x * sin_a + y * cos_a
            femm.mi_addblocklabel(rx, ry)
            femm.mi_selectlabel(rx, ry)
            femm.mi_setblockprop(mat, 0, mesh_size, "<None>", prop_angle, 0, 0)
            femm.mi_clearselected()
            
        place_label(lbl_shaft_air[0], 0, "Air", mesh_shaft)
        place_label(lbl_inner_core[0], 0, rotor_iron_mat, mesh_rotor_core)
        place_label(lbl_central_air[0], 0, "Air", mesh_rotor_core) 
        place_label(lbl_pole_core[0], 0, rotor_iron_mat, mesh_rotor_embedded_pole) 
        place_label(lbl_rib_core[0], 0, rotor_iron_mat, mesh_rib)  
        place_label(ax_out, ay_out, "Air", mesh_rotor_core); place_label(ax_out, -ay_out, "Air", mesh_rotor_core)
        
        polarity_offset = 0 if i % 2 != 0 else 180 
        place_label(xc, yc, magnet_mat, mesh_magnet, mag_angle_base + math.degrees(rot_angle) + polarity_offset)
        place_label(xc, -yc, magnet_mat, mesh_magnet, -mag_angle_base + math.degrees(rot_angle) + polarity_offset)

    femm.mi_zoomnatural()
    femm.mi_saveas(filename)
    femm.mi_close()
    femm.closefemm()
    return 1