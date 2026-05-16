import cmath 
import numpy as np
import femm
import math
import os
import csv
import time
import random
import concurrent.futures
import pandas as pd

# === INITIAL MOTOR PARAMETERS ===
Pm =  10.0  # Mechanical power kW
Vl = math.sqrt(3)*230  # Line-to-Line Voltage in Volts
PF = 0.8  # Power Factor
Ns = 10000 #Speed
n_eff = 0.92 #Efficiency

# === MOTOR TOPOLOGY ===
q = 48 #Number of slots
p = 8  #Number of poles
Np = 2 #Parallel paths
m = 3 #Number of phases
g = 0.8 #Air-gap

# == WINDING PARAMETERS ==
kw = 0.975 #Winding Factor
connection  = "STAR" #Winding connection
winding_layers = "DOUBLE" #Winding layers
gamma_emf = 0.91 #Eph/Vph

# == SLOT PARAMETERS ==
bs0 = 0.2 #Slot opening width
hs0 = 0.3 #Slot opening depth
hs1 = 0.5 #Slot opening dd

B_avg = 2.0/math.pi * 0.9  #Specific Magnetic Loading
ac = 30  #Specific Electric Loading
Jsw = 5 #A/mm2
ar_m = 0.888 #L/D

Bst = 1.7 #Maximum tooth flux density
Bsy = 1.4 #Maximum stator yoke flux density
Bry = 1.4 #Maximum rotor yoke flux density
ki = 0.95 #Stacking factor
kf = 0.35 #Fill factor

# == MECHANICAL & MAGNET CONSTRAINTS (NEW) ==
pole_arc_elec_deg = 120.0 # Electrical pole arc
magnet_width = 12.0 # mm
magnet_thickness = 3.0 # mm

# == INPUTS FEMM ==
BM = np.array([q, p])/math.gcd(q, p) 
winding = np.array([
            [1, 1, -3, -3, 2, 2],
            [1, -3, -3, 2, 2, -1]
          ])

ar = ar_m*p/math.pi 
n_layers = 1 if winding_layers == "SINGLE" else 2 if winding_layers == "DOUBLE" else 0 
nrps = (Ns/60)
G = 1.11*(math.pi**2)*kw*B_avg*ac 

Pe = Pm/n_eff 
P_g = Pe/PF 
f = p*Ns/120 
Il = 1000*Pe/(math.sqrt(3)*Vl*PF)  
Vp = (Vl/math.sqrt(3)) if connection == "STAR" else Vl if connection == "DELTA" else 0 
Ip = Il if connection == "STAR" else (Il/math.sqrt(3)) if connection == "DELTA" else 0 
Nc = (q/m/2) if winding_layers == "SINGLE" else (q/m) if winding_layers == "DOUBLE" else 0 
Nc_p = Nc/Np 
Ic = Ip/Np 
Vc = Vp/Nc_p 

D2L = P_g/G/nrps*1000000000 

#Derived Dimension
D = math.ceil((D2L/ar_m)**(1/3))
L = math.ceil(D*ar_m)
Dr = D - g*2 #Rotor Outer Diameter

phi = B_avg*math.pi*D*L/1000000*1000 
phi_p = phi/p 
phi_st_max = phi_p/2*math.sin(math.pi*p/q) 

wt = phi_st_max/(Bst*L*ki)*1000 
bs1 = 2*(((math.tan(math.pi/q))*(D/2+hs0+hs1))-((wt/2)/(math.cos(math.pi/q))))
wsy = (phi_p/2)/(L/1000*ki*Bsy)

Ep = gamma_emf*Vp 
Nt_ph = Ep/(4.44*f*kw*phi_p)*1000 
Nt_c = Nt_ph/Nc_p 

cAsc = Ic/Jsw 
cAca = cAsc*Nt_c 
gAca = cAca/kf 
sA = n_layers*gAca 

d = (-bs1 + math.sqrt((bs1**2) + (4*math.tan(math.pi/q)*sA)))/(2*math.tan(math.pi/q))
bs2 = bs1 + 2*d*(math.tan(math.pi/q))
Ds = D + 2*(hs0+hs1+d+wsy)


# === HELPER FUNCTIONS ===
def calculate_carters_coefficient(D, q, bs0, g):
    tau_s = (math.pi * D) / q  
    ratio = bs0 / (2 * g)
    gamma = (4 / math.pi) * (ratio * math.atan(ratio) - math.log(math.sqrt(1 + ratio**2)))
    Kc = tau_s / (tau_s - gamma * g)
    return Kc

def get_circle_line_intersection(R, m, c):
    discriminant = R**2 * (1 + m**2) - c**2
    if discriminant < 0:
        raise ValueError("Line does not intersect the circle.")
    x0 = (-m*c + math.sqrt(discriminant)) / (1 + m**2)
    y0 = m * x0 + c
    return x0, y0

def detect_mtpa_peak(gamma_torque_list):
    """Finds the maximum average torque and its corresponding gamma from a list of tuples."""
    best_gamma = 0
    max_torque = -float('inf')
    for g_val, t_val in gamma_torque_list:
        if t_val > max_torque:
            max_torque = t_val
            best_gamma = g_val
    return best_gamma, max_torque


# === GEOMETRY & MESHING FUNCTIONS ===
def drawStator(D, Ds, q, hs0, hs1, d, bs0, bs1, bs2, wsy, BM, winding, Nt_c, filename, 
               mesh_agap=0.15, mesh_stator=None, mesh_coil=None):
    femm.openfemm()
    femm.opendocument(filename)
    
    actual_mesh_stator = mesh_stator if mesh_stator is not None else (wsy / 5)
    actual_mesh_coil = mesh_coil if mesh_coil is not None else min(d / 2, 2.0)
    
    bore_radius = D / 2
    stator_outer_radius = Ds/2
    slot_angle = 360 / q    

    angle_slot_opening = math.asin(bs0 / (2 * bore_radius))
    
    temp0 = complex(bore_radius, 0)* (cmath.exp(1j * angle_slot_opening))
    temp1 = complex(bore_radius, 0)* (cmath.exp(1j * slot_angle/2 * math.pi/180))
    temp2 = complex(stator_outer_radius, 0)* (cmath.exp(1j * slot_angle/2 * math.pi/180))
    
    outer = (temp0.real,temp0.imag)
    top = (temp0.real + hs0 , 0)
    mid = (temp0.real + hs0 + hs1, 0)
    bottom1 = (temp0.real + hs0 + hs1 + d/2, 0)
    bottom2 = (temp0.real + hs0 + hs1 + d/2 + d/2, 0)
    
    Rlayer1 = temp0.real + hs0 + hs1 + d*0.25
    Rlayer2 = temp0.real + hs0 + hs1 + d*0.75
    
    femm.mi_addnode(outer[0], bs0 / 2); femm.mi_addnode(top[0], bs0 / 2)
    femm.mi_addnode(mid[0], bs1 / 2); femm.mi_addnode(bottom1[0], bs1/2 + (bs2-bs1)/2/2)
    femm.mi_addnode(bottom2[0], bs1/2 + 2*(bs2-bs1)/2/2); femm.mi_addnode(temp1.real, temp1.imag)
    femm.mi_addnode(temp2.real, temp2.imag)
    
    femm.mi_addnode(outer[0], -bs0 / 2); femm.mi_addnode(top[0], -bs0 / 2)
    femm.mi_addnode(mid[0], -bs1 / 2); femm.mi_addnode(bottom1[0], -(bs1/2 + (bs2-bs1)/2/2))
    femm.mi_addnode(bottom2[0], -(bs1/2 + 2*(bs2-bs1)/2/2)); femm.mi_addnode(temp1.real, -temp1.imag)
    femm.mi_addnode(temp2.real, -temp2.imag)
    
    femm.mi_addsegment(outer[0], bs0 / 2, top[0], bs0 / 2); femm.mi_addsegment(outer[0], -bs0 / 2, top[0], -bs0 / 2)
    femm.mi_addsegment(top[0], bs0 / 2, mid[0], bs1 / 2); femm.mi_addsegment(top[0], -bs0 / 2, mid[0], -bs1 / 2)
    femm.mi_addsegment(mid[0], -bs1 / 2, mid[0], bs1 / 2)
    femm.mi_addsegment(mid[0], bs1 / 2, bottom1[0], (bs1/2 + (bs2-bs1)/2/2))
    femm.mi_addsegment(mid[0], -bs1 / 2, bottom1[0], -(bs1/2 + (bs2-bs1)/2/2))
    femm.mi_addsegment(bottom1[0], (bs1/2 + (bs2-bs1)/2/2), bottom1[0], -(bs1/2 + (bs2-bs1)/2/2))
    femm.mi_addsegment(bottom1[0], (bs1/2 + (bs2-bs1)/2/2), bottom2[0], bs2 / 2)
    femm.mi_addsegment(bottom1[0], -(bs1/2 + (bs2-bs1)/2/2), bottom2[0], -bs2 / 2)
    femm.mi_addsegment(bottom2[0], bs2 / 2, bottom2[0], -bs2 / 2)
    femm.mi_addsegment(temp1.real, temp1.imag, temp2.real, temp2.imag)
    femm.mi_addsegment(temp1.real, -temp1.imag, temp2.real, -temp2.imag)
    
    femm.mi_addarc(outer[0], bs0 / 2, temp1.real, temp1.imag, slot_angle/2 - angle_slot_opening*2*180/2/math.pi, 1)
    femm.mi_addarc(temp1.real, -temp1.imag, outer[0], -bs0 / 2, slot_angle/2 - angle_slot_opening*2*180/2/math.pi, 1)
    femm.mi_addarc(temp2.real, -temp1.imag, temp2.real, temp1.imag, slot_angle, 1)
    
    femm.mi_selectcircle(0,0,stator_outer_radius*1.5,4)
    femm.mi_setgroup(1001)
    femm.mi_selectgroup(1001)
    femm.mi_moverotate(0,0,slot_angle/2)
    femm.mi_clearselected()
    femm.mi_selectgroup(1001)
    femm.mi_copyrotate(0,0,slot_angle, int(BM[0])-1)
    
    femm.mi_addnode((bore_radius - g/3), 0)
    femm.mi_addsegment(bore_radius, 0, (bore_radius - g/3),0)
    
    femm.mi_addnode((bore_radius - g/3) * np.cos(np.radians(BM[0]*slot_angle)),
                    (bore_radius - g/3) * np.sin(np.radians(BM[0]*slot_angle)))
    femm.mi_addsegment((bore_radius) * np.cos(np.radians(BM[0]*slot_angle)),
                    (bore_radius) * np.sin(np.radians(BM[0]*slot_angle)),
                    (bore_radius - g/3) * np.cos(np.radians(BM[0]*slot_angle)),
                    (bore_radius - g/3) * np.sin(np.radians(BM[0]*slot_angle)))
    
    femm.mi_addarc(bore_radius - g/3, 0, 
                   (bore_radius - g/3) * np.cos(np.radians(BM[0]*slot_angle)),
                   (bore_radius - g/3) * np.sin(np.radians(BM[0]*slot_angle)), BM[0]*slot_angle, 1)
        
    slots_per_pole = len(winding[0])
    for i in range(int(BM[0])):
        yokeLabel = (stator_outer_radius - wsy/2)*np.exp(1j*np.radians(slot_angle/2 + slot_angle*i))
        femm.mi_addblocklabel(yokeLabel.real, yokeLabel.imag)
        femm.mi_selectlabel(yokeLabel.real, yokeLabel.imag)
        femm.mi_setblockprop("US Steel Type 2-S 0.024 inch thickness", 0, actual_mesh_stator)
        femm.mi_clearselected()
        
        w_idx = i % slots_per_pole
        pole_idx = i // slots_per_pole
        polarity = 1 if pole_idx % 2 == 0 else -1
        w1 = winding[0, w_idx] * polarity
        w2 = winding[1, w_idx] * polarity
        
        Lay1 = Rlayer1*np.exp(1j*np.radians(slot_angle/2 + slot_angle*i))
        femm.mi_addblocklabel(Lay1.real, Lay1.imag)
        femm.mi_selectlabel(Lay1.real, Lay1.imag)
        femm.mi_setblockprop("24 SWG", 0, actual_mesh_coil, f"fase{int(np.abs(w1))}", 0, 1001, Nt_c*(w1/np.abs(w1)))
        femm.mi_clearselected()
        
        Lay2 = Rlayer2*np.exp(1j*np.radians(slot_angle/2 + slot_angle*i))
        femm.mi_addblocklabel(Lay2.real, Lay2.imag)
        femm.mi_selectlabel(Lay2.real, Lay2.imag)
        femm.mi_setblockprop("24 SWG", 0, actual_mesh_coil, f"fase{int(np.abs(w2))}", 0, 1001, Nt_c*(w2/np.abs(w2)))
        femm.mi_clearselected()
        
    femm.mi_selectcircle(0,0,stator_outer_radius*1.5,4)
    femm.mi_setgroup(1001)
    femm.mi_clearselected()
    
    stator_air_radius = bore_radius - (g / 6)
    stator_air_angle = math.radians(slot_angle / 2)
    s_air_x = stator_air_radius * math.cos(stator_air_angle)
    s_air_y = stator_air_radius * math.sin(stator_air_angle)
    femm.mi_addblocklabel(s_air_x, s_air_y)
    femm.mi_selectlabel(s_air_x, s_air_y)
    femm.mi_setblockprop("Air", 0, mesh_agap, "<None>", 0, 0, 0)
    femm.mi_clearselected()
    
    femm.mi_saveas(filename)
    femm.mi_close()
    femm.closefemm()
    return 1

def drawFullRotorIPM(Dr, tr_t, g, p, BM, pole_arc_elec_deg, magnet_thickness, magnet_width, rotor_iron_mat, magnet_mat, filename, 
                     mesh_agap=0.15, mesh_rotor_core=1.0, mesh_rotor_embedded_pole=0.75, mesh_rib=0.15, mesh_magnet=0.5, mesh_shaft=5.0):
    femm.openfemm()
    femm.opendocument(filename)
    
    rotor_radius = Dr / 2
    shaft_radius = rotor_radius * 0.2 
    r_air = rotor_radius + (g / 3)
    
    theta_pole = 360 / p
    ang_pole_rad = math.radians(theta_pole / 2)
    pole_arc_mech_deg = pole_arc_elec_deg / (p / 2)
    theta_arc_rad = math.radians(pole_arc_mech_deg / 2)
    
    fixed_tr_t = 2.0
    ref_design_radius = rotor_radius - fixed_tr_t
    
    x1_ref = ref_design_radius * math.cos(theta_arc_rad)
    y1_ref = ref_design_radius * math.sin(theta_arc_rad)
    
    pocket_length = 1.2 * magnet_width
    if pocket_length <= y1_ref: pocket_length = y1_ref * 1.01 
        
    beta_rad = math.asin(y1_ref / pocket_length)
    m1 = math.tan(beta_rad) 
    
    x2 = x1_ref - (y1_ref / m1)
    shift_distance = magnet_thickness / math.sin(beta_rad)
    x3 = x2 - shift_distance
    
    x_mid_in = (x2 + x3) / 2; y_mid_in = 0
    c_mid = -m1 * x_mid_in
    x_mid_out, y_mid_out = get_circle_line_intersection(ref_design_radius, m1, c_mid)
    
    xc = (x_mid_in + x_mid_out) / 2; yc = (y_mid_out) / 2
    dx = x_mid_out - x_mid_in; dy = y_mid_out - y_mid_in
    actual_pocket_len = math.hypot(dx, dy)
    
    ux = dx / actual_pocket_len; uy = dy / actual_pocket_len
    nx = -uy; ny = ux   
    if (nx*xc + ny*yc) < 0:
        nx = uy; ny = -ux
        
    hl = magnet_width / 2; ht = magnet_thickness / 2
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
    
    femm.mi_addnode(x1, y1); femm.mi_addnode(x4, y4)
    femm.mi_addnode(x1, -y1); femm.mi_addnode(x4, -y4)
    femm.mi_addnode(x2, 0); femm.mi_addnode(x3, 0)
    
    femm.mi_addsegment(0, 0, x_sh_p, y_sh_p); femm.mi_addsegment(x_sh_p, y_sh_p, x_des_p, y_des_p)
    femm.mi_addsegment(x_des_p, y_des_p, x_out_p, y_out_p)
    femm.mi_addsegment(0, 0, x_sh_p, -y_sh_p); femm.mi_addsegment(x_sh_p, -y_sh_p, x_des_p, -y_des_p)
    femm.mi_addsegment(x_des_p, -y_des_p, x_out_p, -y_out_p)

    for flip_y in [1, -1]:
        mx1_f = mx1; my1_f = my1 * flip_y
        mx2_f = mx2; my2_f = my2 * flip_y
        mx3_f = mx3; my3_f = my3 * flip_y
        mx4_f = mx4; my4_f = my4 * flip_y
        femm.mi_addnode(mx1_f, my1_f); femm.mi_addnode(mx2_f, my2_f)
        femm.mi_addnode(mx3_f, my3_f); femm.mi_addnode(mx4_f, my4_f)
        femm.mi_addsegment(mx1_f, my1_f, mx2_f, my2_f); femm.mi_addsegment(mx2_f, my2_f, mx3_f, my3_f)
        femm.mi_addsegment(mx3_f, my3_f, mx4_f, my4_f); femm.mi_addsegment(mx4_f, my4_f, mx1_f, my1_f)

    femm.mi_addsegment(x1, y1, x2, 0); femm.mi_addsegment(x4, y4, x3, 0)
    femm.mi_addsegment(x1, -y1, x2, 0); femm.mi_addsegment(x4, -y4, x3, 0)

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

    total_angle_deg = int(BM[1]) * theta_pole
    n_arc_chunks = math.ceil(total_angle_deg / 179.0) 
    chunk_deg = total_angle_deg / n_arc_chunks
    
    for i in range(n_arc_chunks):
        start_rad = math.radians(i * chunk_deg)
        end_rad = math.radians((i + 1) * chunk_deg)
        xa = r_air * math.cos(start_rad); ya = r_air * math.sin(start_rad)
        xb = r_air * math.cos(end_rad); yb = r_air * math.sin(end_rad)
        femm.mi_addnode(xa, ya); femm.mi_addnode(xb, yb)
        femm.mi_addarc(xa, ya, xb, yb, chunk_deg, 1)

    if int(BM[1]) < p:
        total_rad = math.radians(total_angle_deg)
        xs_out = rotor_radius; ys_out = 0; xs_air = r_air; ys_air = 0
        femm.mi_addnode(xs_out, ys_out); femm.mi_addnode(xs_air, ys_air)
        femm.mi_addsegment(xs_out, ys_out, xs_air, ys_air)
        xe_out = rotor_radius * math.cos(total_rad); ye_out = rotor_radius * math.sin(total_rad)
        xe_air = r_air * math.cos(total_rad); ye_air = r_air * math.sin(total_rad)
        femm.mi_addnode(xe_out, ye_out); femm.mi_addnode(xe_air, ye_air)
        femm.mi_addsegment(xe_out, ye_out, xe_air, ye_air)

    mag_angle_base = math.degrees(math.atan2(ny, nx))
    lbl_shaft_air = (shaft_radius / 2, 0)
    lbl_inner_core = ((shaft_radius + x3) / 2, 0)
    lbl_central_air = ((x3 + x2) / 2, 0)
    lbl_pole_core = ((x2 + design_radius) / 2, 0)
    lbl_rib_core = ((design_radius + rotor_radius) / 2, 0)
    
    ax_out = (x_mid_out + (xc + hl*ux)) / 2
    ay_out = (y_mid_out + (yc + hl*uy)) / 2
    
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
            rx = x * cos_a - y * sin_a
            ry = x * sin_a + y * cos_a
            femm.mi_addblocklabel(rx, ry)
            femm.mi_selectlabel(rx, ry)
            femm.mi_setblockprop(mat, 0, mesh_size, "<None>", prop_angle, 0, 0)
            femm.mi_clearselected()
            
        place_label(lbl_shaft_air[0], 0, "Air", mesh_shaft)
        place_label(lbl_inner_core[0], 0, rotor_iron_mat, mesh_rotor_core)
        place_label(lbl_central_air[0], 0, "Air", mesh_rotor_core) 
        place_label(lbl_pole_core[0], 0, rotor_iron_mat, mesh_rotor_embedded_pole) 
        place_label(lbl_rib_core[0], 0, rotor_iron_mat, mesh_rib)  
        
        place_label(ax_out, ay_out, "Air", mesh_rotor_core)
        place_label(ax_out, -ay_out, "Air", mesh_rotor_core)
        
        polarity_offset = 180 if i % 2 != 0 else 0 
        
        mag_angle_top = mag_angle_base + math.degrees(rot_angle) + polarity_offset
        place_label(xc, yc, magnet_mat, mesh_magnet, mag_angle_top)
        mag_angle_bot = -mag_angle_base + math.degrees(rot_angle) + polarity_offset
        place_label(xc, -yc, magnet_mat, mesh_magnet, mag_angle_bot)

    femm.mi_zoomnatural()
    femm.mi_saveas(filename)
    femm.mi_close()
    femm.closefemm()
    return 1

def assign_all_boundaries(D, Ds, Dr, tr_t, g, p, q, BM, filename):
    femm.openfemm()
    femm.opendocument(filename)
    
    is_anti_periodic = (int(BM[1]) % 2 != 0)
    b_type = 5 if is_anti_periodic else 4 
    agap_type = 7 if is_anti_periodic else 6 
    
    theta_pole = 360 / p
    slot_angle = 360 / q
    inner_angle_deg = int(BM[1]) * theta_pole  
    outer_angle_deg = int(BM[0]) * slot_angle  
    
    femm.mi_addboundprop("A=0", 0, 0, 0, 0, 0, 0, 0, 0, 0)
    femm.mi_addboundprop("AGap", 0, 0, 0, 0, 0, 0, inner_angle_deg, outer_angle_deg, agap_type) 
    
    for i in range(1, 5):
        femm.mi_addboundprop(f"r{i}", 0, 0, 0, 0, 0, 0, 0, 0, b_type)
    for i in range(1, 3):
        femm.mi_addboundprop(f"s{i}", 0, 0, 0, 0, 0, 0, 0, 0, b_type)

    rotor_radius = Dr / 2
    shaft_radius = rotor_radius * 0.2
    design_radius = rotor_radius - tr_t
    r_air = rotor_radius + (g / 3)
    bore_radius = D / 2
    stator_outer_radius = Ds / 2
    s_air = bore_radius - (g / 3)
    
    total_r_rad = math.radians(inner_angle_deg)
    total_s_rad = math.radians(outer_angle_deg)
    
    r_mids = [shaft_radius / 2, (shaft_radius + design_radius) / 2, (design_radius + rotor_radius) / 2, (rotor_radius + r_air) / 2]
    for i, r_mid in enumerate(r_mids):
        b_name = f"r{i+1}"
        femm.mi_selectsegment(r_mid, 0)
        femm.mi_setsegmentprop(b_name, 0, 1, 0, 0)
        femm.mi_clearselected()
        femm.mi_selectsegment(r_mid * math.cos(total_r_rad), r_mid * math.sin(total_r_rad))
        femm.mi_setsegmentprop(b_name, 0, 1, 0, 0)
        femm.mi_clearselected()

    s_mids = [(bore_radius + s_air) / 2, (bore_radius + stator_outer_radius) / 2]
    for i, s_mid in enumerate(s_mids):
        b_name = f"s{i+1}"
        femm.mi_selectsegment(s_mid, 0)
        femm.mi_setsegmentprop(b_name, 0, 1, 0, 0)
        femm.mi_clearselected()
        femm.mi_selectsegment(s_mid * math.cos(total_s_rad), s_mid * math.sin(total_s_rad))
        femm.mi_setsegmentprop(b_name, 0, 1, 0, 0)
        femm.mi_clearselected()

    for i in range(int(BM[0])):
        theta_arc = math.radians((i + 0.5) * slot_angle)
        x_arc = stator_outer_radius * math.cos(theta_arc)
        y_arc = stator_outer_radius * math.sin(theta_arc)
        femm.mi_selectarcsegment(x_arc, y_arc)
        femm.mi_setarcsegmentprop(1, "A=0", 0, 0)
        femm.mi_clearselected()

    femm.mi_selectarcsegment(s_air * math.cos(total_s_rad / 2), s_air * math.sin(total_s_rad / 2))
    femm.mi_setarcsegmentprop(1, "AGap", 0, 0)
    femm.mi_clearselected()

    n_arc_chunks = math.ceil(inner_angle_deg / 179.0) 
    chunk_deg = inner_angle_deg / n_arc_chunks
    for i in range(n_arc_chunks):
        start_rad = math.radians(i * chunk_deg)
        mid_rad = start_rad + math.radians(chunk_deg) / 2
        arc_mid_x = r_air * math.cos(mid_rad)
        arc_mid_y = r_air * math.sin(mid_rad)
        femm.mi_selectarcsegment(arc_mid_x, arc_mid_y)
        femm.mi_setarcsegmentprop(1, "AGap", 0, 0)
        femm.mi_clearselected()

    femm.mi_zoomnatural()
    femm.mi_saveas(filename)
    femm.mi_close()
    femm.closefemm()
    print("All cyclic, outer, and AGap boundaries successfully assigned!")
    return 1


# === SIMULATION & MTPA FUNCTIONS ===
def parallel_analyze_step(args):
    step_num, theta_mech_deg, i_a, i_b, i_c, p, BM, base_filename, step_filename, band_name, save_files = args
    time.sleep(random.uniform(0.1, 1.5))
    
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
        
        sym_mult = p / BM[1]
        raw_torque = femm.mo_gapintegral(band_name, 0)
        torque_total = raw_torque * sym_mult
        
        props_a = femm.mo_getcircuitproperties("fase1")
        props_b = femm.mo_getcircuitproperties("fase2")
        props_c = femm.mo_getcircuitproperties("fase3")
        
        v_a = props_a[1] * sym_mult; flux_a = props_a[2] * sym_mult
        v_b = props_b[1] * sym_mult; flux_b = props_b[2] * sym_mult
        v_c = props_c[1] * sym_mult; flux_c = props_c[2] * sym_mult
        
        femm.mo_close()
        femm.mi_close()
        femm.closefemm()
        
        # Clean up heavy FEMM files instantly if flag is off
        if not save_files:
            try:
                if os.path.exists(step_filename): os.remove(step_filename)
                ans_file = step_filename.replace('.fem', '.ans')
                if os.path.exists(ans_file): os.remove(ans_file)
            except:
                pass
                
        return {
            "step_index": step_num, "Torque_Nm": torque_total,
            "I_A": props_a[0], "I_B": props_b[0], "I_C": props_c[0],
            "V_drop_A": v_a, "V_drop_B": v_b, "V_drop_C": v_c,
            "Flux_A_Wb": flux_a, "Flux_B_Wb": flux_b, "Flux_C_Wb": flux_c
        }
    except Exception as e:
        print(f"ERROR on step {step_num}: {e}")
        try: femm.closefemm()
        except: pass
        return None
    
def run_motor_simulation(I_rms, gamma_elec_deg, initial_pos, Ns_rpm, p, BM, theta_start_elec, theta_end_elec, num_steps, filename, base_folder, band_name="AGap", save_files=True):
    sim_folder_name = f"Sim_{I_rms:.2f}A_{gamma_elec_deg}deg_{num_steps}steps"
    full_sim_folder = os.path.join(base_folder, sim_folder_name)
    steps_folder = os.path.join(full_sim_folder, "steps") if save_files else full_sim_folder
    
    if not os.path.exists(steps_folder):
        os.makedirs(steps_folder) 
        
    results = {
        "Theta_Elec_deg": [], "Theta_Mech_deg": [], "Torque_Nm": [],
        "I_A": [], "I_B": [], "I_C": [],
        "V_drop_A": [], "V_drop_B": [], "V_drop_C": [],
        "Flux_A_Wb": [], "Flux_B_Wb": [], "Flux_C_Wb": []
    }
    
    I_peak = I_rms * math.sqrt(2)
    step_size_elec = (theta_end_elec - theta_start_elec) / (num_steps - 1) if num_steps > 1 else 0
    
    tasks = []
    for step in range(num_steps):
        theta_elec = theta_start_elec + (step * step_size_elec)
        theta_mech = theta_elec / (p / 2) + initial_pos
        
        i_a = I_peak * math.cos(math.radians(theta_elec + gamma_elec_deg))
        i_b = I_peak * math.cos(math.radians(theta_elec - 120 + gamma_elec_deg))
        i_c = I_peak * math.cos(math.radians(theta_elec + 120 + gamma_elec_deg))
        
        step_filename = os.path.join(steps_folder, f"step_{step:03d}_{gamma_elec_deg}deg.fem")
        tasks.append((step, theta_mech, i_a, i_b, i_c, p, BM, filename, step_filename, band_name, save_files))
        
        results["Theta_Elec_deg"].append(theta_elec)
        results["Theta_Mech_deg"].append(theta_mech)

    max_cores = max(1, min(6, os.cpu_count() - 1)) 
    with concurrent.futures.ProcessPoolExecutor(max_workers=max_cores) as executor:
        unordered_results = list(executor.map(parallel_analyze_step, tasks))

    valid_results = [r for r in unordered_results if r is not None]
    sorted_results = sorted(valid_results, key=lambda x: x["step_index"])
    
    for r in sorted_results:
        results["Torque_Nm"].append(r["Torque_Nm"])
        results["I_A"].append(r["I_A"]); results["I_B"].append(r["I_B"]); results["I_C"].append(r["I_C"])
        results["V_drop_A"].append(r["V_drop_A"]); results["V_drop_B"].append(r["V_drop_B"]); results["V_drop_C"].append(r["V_drop_C"])
        results["Flux_A_Wb"].append(r["Flux_A_Wb"]); results["Flux_B_Wb"].append(r["Flux_B_Wb"]); results["Flux_C_Wb"].append(r["Flux_C_Wb"])

    if save_files:
        csv_file_path = os.path.join(full_sim_folder, "simulation_results.csv")
        keys = results.keys()
        with open(csv_file_path, 'w', newline='') as output_file:
            dict_writer = csv.writer(output_file)
            dict_writer.writerow(keys)
            dict_writer.writerows(zip(*[results[key] for key in keys]))
            
    avg_torque = np.mean(results["Torque_Nm"])
    return results, avg_torque
# === THE TWO-STAGE MTPA ALGORITHM ===
def run_mtpa_sweep(I_min, I_max, I_steps, gamma_min, gamma_max, gamma_steps, enable_refinement, initial_pos, Ns_rpm, p, BM, theta_start, theta_end, num_steps, base_filename, rib_folder, **kwargs):
    
    # 1. Input Validation for Refinement
    if enable_refinement:
        if 'fine_window' not in kwargs or 'fine_step' not in kwargs:
            raise ValueError("ERROR: 'enable_refinement' is True, but 'fine_window' and/or 'fine_step' were not provided in kwargs.")
        fine_window = kwargs['fine_window']
        fine_step = kwargs['fine_step']
    
    mtpa_all_samples = []
    mtpa_optimal_curve = []
    
    # 2. Generate Base Arrays
    current_array = np.linspace(I_min, I_max, I_steps)
    coarse_gammas = np.linspace(gamma_min, gamma_max, gamma_steps)
    
    for current in current_array:
        if current <= 0.01:
            mtpa_optimal_curve.append({"I_rms_A": 0.0, "Optimal_Gamma_deg": 0.0, "Max_Avg_Torque_Nm": 0.0})
            continue
            
        print(f"\n   [MTPA] Analyzing Current = {current:.2f} A")
        gamma_torque_current = []
        simulated_gammas = set()
        
        # --- STAGE 1: COARSE SWEEP ---
        for g_val in coarse_gammas:
            print(f"      -> Coarse Map: Gamma = {g_val:.1f}°")
            _, avg_t = run_motor_simulation(current, g_val, initial_pos, Ns_rpm, p, BM, theta_start, theta_end, num_steps, base_filename, rib_folder, save_files=False)
            gamma_torque_current.append((g_val, avg_t))
            simulated_gammas.add(g_val)
            mtpa_all_samples.append({"I_rms_A": current, "Gamma_deg": g_val, "Avg_Torque_Nm": avg_t, "Type": "Coarse"})
            
        best_coarse_g, _ = detect_mtpa_peak(gamma_torque_current)
        print(f"      -> Peak detected near {best_coarse_g:.1f}°")
        
        # --- STAGE 2: FINE REFINEMENT ---
        if enable_refinement:
            fine_start = max(gamma_min, best_coarse_g - fine_window)
            fine_end = min(gamma_max, best_coarse_g + fine_window)
            fine_gammas = np.arange(fine_start, fine_end + fine_step, fine_step)
            
            for fg in fine_gammas:
                # Due to floating point precision, check if fg is very close to an already simulated gamma
                already_simulated = any(math.isclose(fg, sg, abs_tol=1e-5) for sg in simulated_gammas)
                if already_simulated:
                    continue
                    
                print(f"      -> Fine Map: Gamma = {fg:.1f}°")
                _, avg_t = run_motor_simulation(current, fg, initial_pos, Ns_rpm, p, BM, theta_start, theta_end, num_steps, base_filename, rib_folder, save_files=False)
                gamma_torque_current.append((fg, avg_t))
                simulated_gammas.add(fg)
                mtpa_all_samples.append({"I_rms_A": current, "Gamma_deg": fg, "Avg_Torque_Nm": avg_t, "Type": "Fine"})
                
        # --- FINAL DETECTION ---
        true_best_g, true_max_t = detect_mtpa_peak(gamma_torque_current)
        mtpa_optimal_curve.append({"I_rms_A": current, "Optimal_Gamma_deg": true_best_g, "Max_Avg_Torque_Nm": true_max_t})
        print(f"   [MTPA] Locked Peak for {current:.2f} A -> Gamma: {true_best_g:.1f}°, Torque: {true_max_t:.4f} Nm")
        
    return mtpa_all_samples, mtpa_optimal_curve



# === MASTER EXECUTION BLOCK ===
if __name__ == '__main__':
    rotor_iron_mat = "US Steel Type 2-S 0.024 inch thickness"
    magnet_mat = "N35"

    rib_thicknesses = [0.8]
    master_sweep_folder = "Parametric_Rib_Sweep"

    # Simulation Resolution (60 deg elec in 19 steps = exactly 1 full 6th harmonic ripple cycle)
    mtpa_theta_start = 0.0
    mtpa_theta_end = 60.0 
    mtpa_num_steps = 5  
    
    # Sweep Setup
    I_min, I_max, I_steps = 0.0, Ic, 2
    gamma_min, gamma_max, gamma_steps = 0.0, 90.0, 5
    
    print("\n" + "="*60)
    print("STARTING PARAMETRIC GEOMETRY & MTPA SWEEP")
    print("="*60)

    for current_tr_t in rib_thicknesses:
        rib_folder = os.path.join(master_sweep_folder, f"Rib_{current_tr_t}mm")
        if not os.path.exists(rib_folder): os.makedirs(rib_folder)
        current_filename = os.path.join(rib_folder, f"Machine_Rib_{current_tr_t}mm.fem")

        print(f"\n[{current_tr_t} mm] STEP 1: Building Geometry...")
        
        femm.openfemm()
        femm.newdocument(0)
        femm.mi_probdef(0, "millimeters", "planar", 1e-8, L, 30)
        femm.mi_getmaterial("US Steel Type 2-S 0.024 inch thickness")
        femm.mi_getmaterial("Air")
        femm.mi_getmaterial("24 SWG")
        femm.mi_getmaterial("N35") 
        femm.mi_smartmesh(0)
        for i in range(m): femm.mi_addcircprop(f"fase{i+1}", 0, 1)

        femm.mi_saveas(current_filename)
        femm.mi_close()
        femm.closefemm()

        drawStator(D, Ds, q, hs0, hs1, d, bs0, bs1, bs2, wsy, BM, winding, Nt_c, current_filename, 
                mesh_agap=0.1, mesh_stator=0.75, mesh_coil=1.5)               
        drawFullRotorIPM(Dr, current_tr_t, g, p, BM, pole_arc_elec_deg, magnet_thickness, magnet_width, rotor_iron_mat, magnet_mat, current_filename, 
                        mesh_agap=0.1, mesh_rotor_core=0.75, mesh_rotor_embedded_pole=0.1, mesh_rib=0.10, mesh_magnet=0.3, mesh_shaft=5.0)
        
        assign_all_boundaries(D, Ds, Dr, current_tr_t, g, p, q, BM, current_filename)
        print(f"[{current_tr_t} mm] Geometry generated successfully!")

        print(f"\n[{current_tr_t} mm] STEP 2: Running Adaptive MTPA Algorithm...")
        
        all_samples, optimal_curve = run_mtpa_sweep(
            I_min=I_min, I_max=I_max, I_steps=I_steps,
            gamma_min=gamma_min, gamma_max=gamma_max, gamma_steps=gamma_steps,
            enable_refinement=True,
            fine_window=5, fine_step=2, # Passed as kwargs!
            initial_pos=3.75,
            Ns_rpm=Ns,
            p=p,
            BM=BM,
            theta_start=mtpa_theta_start,
            theta_end=mtpa_theta_end,
            num_steps=mtpa_num_steps,
            base_filename=current_filename,
            rib_folder=rib_folder
        )
        
        # Save Outputs
        pd.DataFrame(all_samples).to_csv(os.path.join(rib_folder, "MTPA_All_Samples.csv"), index=False)
        pd.DataFrame(optimal_curve).to_csv(os.path.join(rib_folder, "MTPA_Optimal_Curve.csv"), index=False)
        print(f"[{current_tr_t} mm] MTPA Sweep Complete & Data Saved!")

    print("\n" + "="*60)
    print("ALL PARAMETRIC SWEEPS COMPLETED SUCCESSFULLY!")
    print("="*60)