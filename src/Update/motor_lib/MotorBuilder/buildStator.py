import cmath 
import femm
import math

def draw_stator(
    *,
    D, Ds, p, q, g, hs0, hs1, d, bs0, bs1, bs2, wsy, BM, winding, Nt_c, filename, 
    stator_iron_mat, coil_mat, mesh_agap, mesh_stator, mesh_coil
):
    femm.openfemm()
    femm.opendocument(filename)
    
    actual_mesh_stator = mesh_stator if mesh_stator is not None else (wsy / 5)
    actual_mesh_coil = mesh_coil if mesh_coil is not None else min(d / 2, 2.0)
    
    bore_radius = D / 2
    stator_outer_radius = Ds / 2
    slot_angle = 360 / q    

    angle_slot_opening = math.asin(bs0 / (2 * bore_radius))
    
    temp0 = complex(bore_radius, 0) * cmath.exp(1j * angle_slot_opening)
    temp1 = complex(bore_radius, 0) * cmath.exp(1j * math.radians(slot_angle / 2))
    temp2 = complex(stator_outer_radius, 0) * cmath.exp(1j * math.radians(slot_angle / 2))
    
    outer_x = temp0.real; top_x = temp0.real + hs0; mid_x = temp0.real + hs0 + hs1
    bottom1_x = temp0.real + hs0 + hs1 + d/2; bottom2_x = temp0.real + hs0 + hs1 + d
    
    Rlayer1 = temp0.real + hs0 + hs1 + (d * 0.25)
    Rlayer2 = temp0.real + hs0 + hs1 + (d * 0.75)
    
    # --- DRAW NODES & SEGMENTS ---
    femm.mi_addnode(outer_x, bs0 / 2); femm.mi_addnode(top_x, bs0 / 2)
    femm.mi_addnode(mid_x, bs1 / 2); femm.mi_addnode(bottom1_x, bs1/2 + (bs2-bs1)/4)
    femm.mi_addnode(bottom2_x, bs2 / 2); femm.mi_addnode(temp1.real, temp1.imag)
    femm.mi_addnode(temp2.real, temp2.imag)
    
    femm.mi_addnode(outer_x, -bs0 / 2); femm.mi_addnode(top_x, -bs0 / 2)
    femm.mi_addnode(mid_x, -bs1 / 2); femm.mi_addnode(bottom1_x, -(bs1/2 + (bs2-bs1)/4))
    femm.mi_addnode(bottom2_x, -bs2 / 2); femm.mi_addnode(temp1.real, -temp1.imag)
    femm.mi_addnode(temp2.real, -temp2.imag)
    
    femm.mi_addsegment(outer_x, bs0/2, top_x, bs0/2); femm.mi_addsegment(outer_x, -bs0/2, top_x, -bs0/2)
    femm.mi_addsegment(top_x, bs0/2, mid_x, bs1/2); femm.mi_addsegment(top_x, -bs0/2, mid_x, -bs1/2)
    femm.mi_addsegment(mid_x, -bs1/2, mid_x, bs1/2)
    femm.mi_addsegment(mid_x, bs1/2, bottom1_x, bs1/2 + (bs2-bs1)/4)
    femm.mi_addsegment(mid_x, -bs1/2, bottom1_x, -(bs1/2 + (bs2-bs1)/4))
    femm.mi_addsegment(bottom1_x, bs1/2 + (bs2-bs1)/4, bottom1_x, -(bs1/2 + (bs2-bs1)/4))
    femm.mi_addsegment(bottom1_x, bs1/2 + (bs2-bs1)/4, bottom2_x, bs2/2)
    femm.mi_addsegment(bottom1_x, -(bs1/2 + (bs2-bs1)/4), bottom2_x, -bs2/2)
    femm.mi_addsegment(bottom2_x, bs2/2, bottom2_x, -bs2/2)
    
    femm.mi_addsegment(temp1.real, temp1.imag, temp2.real, temp2.imag)
    femm.mi_addsegment(temp1.real, -temp1.imag, temp2.real, -temp2.imag)
    
    femm.mi_addarc(outer_x, bs0/2, temp1.real, temp1.imag, slot_angle/2 - math.degrees(angle_slot_opening), 1)
    femm.mi_addarc(temp1.real, -temp1.imag, outer_x, -bs0/2, slot_angle/2 - math.degrees(angle_slot_opening), 1)
    femm.mi_addarc(temp2.real, -temp1.imag, temp2.real, temp1.imag, slot_angle, 1)
    
    # --- COPY / ROTATE SLOT ---
    femm.mi_selectcircle(0, 0, stator_outer_radius * 1.5, 4)
    femm.mi_setgroup(1001)
    femm.mi_selectgroup(1001)
    femm.mi_moverotate(0, 0, slot_angle / 2) 
    femm.mi_clearselected()
    
    num_stator_slots = int(BM[0])
    femm.mi_selectgroup(1001)
    femm.mi_copyrotate(0, 0, slot_angle, num_stator_slots - 1)
    femm.mi_clearselected()
    
    # --- DRAW AIRGAP SLIDING BAND ---
    airgap_boundary = bore_radius - (g / 3)
    total_s_angle_deg = num_stator_slots * slot_angle
    
    if num_stator_slots < q:
        # Standard Sector
        femm.mi_addnode(airgap_boundary, 0)
        ag_x = airgap_boundary * math.cos(math.radians(total_s_angle_deg))
        ag_y = airgap_boundary * math.sin(math.radians(total_s_angle_deg))
        femm.mi_addnode(ag_x, ag_y)
        femm.mi_addarc(airgap_boundary, 0, ag_x, ag_y, total_s_angle_deg, 1)
        
        # Stator Radial Boundaries
        b_x = bore_radius * math.cos(math.radians(total_s_angle_deg))
        b_y = bore_radius * math.sin(math.radians(total_s_angle_deg))
        femm.mi_addsegment(bore_radius, 0, airgap_boundary, 0)
        femm.mi_addsegment(b_x, b_y, ag_x, ag_y)
    else:
        # Full Machine Fallback
        femm.mi_addnode(airgap_boundary, 0); femm.mi_addnode(-airgap_boundary, 0)
        femm.mi_addarc(airgap_boundary, 0, -airgap_boundary, 0, 180, 1)
        femm.mi_addarc(-airgap_boundary, 0, airgap_boundary, 0, 180, 1)
    
    # --- BOUNDARIES (PERIODIC & OUTER) ---
    is_anti_periodic = (int(BM[1]) % 2 != 0)
    b_type = 5 if is_anti_periodic else 4 
    
    try:
        femm.mi_addboundprop("A=0", 0, 0, 0, 0, 0, 0, 0, 0, 0)
        if num_stator_slots < q:
            for i in range(1, 3): femm.mi_addboundprop(f"s{i}", 0, 0, 0, 0, 0, 0, 0, 0, b_type)
    except:
        pass

    for i in range(num_stator_slots):
        arc_x = stator_outer_radius * math.cos(math.radians((i + 0.5) * slot_angle))
        arc_y = stator_outer_radius * math.sin(math.radians((i + 0.5) * slot_angle))
        femm.mi_selectarcsegment(arc_x, arc_y)
        femm.mi_setarcsegmentprop(1, "A=0", 0, 0)
        femm.mi_clearselected()

    if num_stator_slots < q:
        s_mids = [(bore_radius + airgap_boundary) / 2, (bore_radius + stator_outer_radius) / 2]
        for i, s_mid in enumerate(s_mids):
            b_name = f"s{i+1}"
            femm.mi_selectsegment(s_mid, 0)
            femm.mi_setsegmentprop(b_name, 0, 1, 0, 0)
            femm.mi_clearselected()
            
            x_end = s_mid * math.cos(math.radians(total_s_angle_deg))
            y_end = s_mid * math.sin(math.radians(total_s_angle_deg))
            femm.mi_selectsegment(x_end, y_end)
            femm.mi_setsegmentprop(b_name, 0, 1, 0, 0)
            femm.mi_clearselected()
    else:
        print("WARNING (Stator): Full machine modeled. Cyclic boundaries skipped.")

    # --- MATERIAL & BLOCK LABEL ASSIGNMENT ---
    slots_per_pole = len(winding[0])
    for i in range(num_stator_slots):
        yokeLabel = (stator_outer_radius - wsy/2) * cmath.exp(1j * math.radians(slot_angle/2 + slot_angle*i))
        femm.mi_addblocklabel(yokeLabel.real, yokeLabel.imag)
        femm.mi_selectlabel(yokeLabel.real, yokeLabel.imag)
        femm.mi_setblockprop(stator_iron_mat, 0, actual_mesh_stator, "<None>", 0, 0, 0)
        femm.mi_clearselected()
        
        w_idx = i % slots_per_pole
        pole_idx = i // slots_per_pole
        polarity = 1 if pole_idx % 2 == 0 else -1
        w1 = winding[0, w_idx] * polarity; w2 = winding[1, w_idx] * polarity
        
        Lay1 = Rlayer1 * cmath.exp(1j * math.radians(slot_angle/2 + slot_angle*i))
        femm.mi_addblocklabel(Lay1.real, Lay1.imag)
        femm.mi_selectlabel(Lay1.real, Lay1.imag)
        femm.mi_setblockprop(coil_mat, 0, actual_mesh_coil, f"fase{int(abs(w1))}", 0, 1001, Nt_c * (w1 / abs(w1)))
        femm.mi_clearselected()
        
        Lay2 = Rlayer2 * cmath.exp(1j * math.radians(slot_angle/2 + slot_angle*i))
        femm.mi_addblocklabel(Lay2.real, Lay2.imag)
        femm.mi_selectlabel(Lay2.real, Lay2.imag)
        femm.mi_setblockprop(coil_mat, 0, actual_mesh_coil, f"fase{int(abs(w2))}", 0, 1001, Nt_c * (w2 / abs(w2)))
        femm.mi_clearselected()
        
    s_air_x = (bore_radius - (g / 6)) * math.cos(math.radians(slot_angle / 2))
    s_air_y = (bore_radius - (g / 6)) * math.sin(math.radians(slot_angle / 2))
    femm.mi_addblocklabel(s_air_x, s_air_y)
    femm.mi_selectlabel(s_air_x, s_air_y)
    femm.mi_setblockprop("Air", 0, mesh_agap, "<None>", 0, 0, 0)
    femm.mi_clearselected()
    
    femm.mi_zoomnatural()
    femm.mi_saveas(filename)
    femm.mi_close()
    femm.closefemm()
    return 1