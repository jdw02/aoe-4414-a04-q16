# ecef_to_sez.py

# Usage: py ecef_to_sez.py o_x_km o_y_km o_z_km x_km y_km z_km
#      - o_ represents the ECEF origin of the SEZ frame(station/observer), 
#        and the other coordinates represent the ECEF position(object position)   

# Parameters:
#  arg1: o_x_km <--- ECEF x origin of SEZ frame (observer)
#  arg2: o_y_km <--- ECEF y origin of SEZ frame (observer)
#  arg3: o_z_km <--- ECEF z origin of SEZ frame (observer)
#  arg4: x_km <--- ECEF x object position 
#  arg5: y_km <--- ECEF y object position 
#  arg6: z_km <--- ECEF z object position 

# Output:
#  Takes ECEF coordinates and outputs SEZ coordinates

# Written by Jayden Warren
# Other contributors: None

# import Python modules
import sys # argv
import math

# "constants"
R_E_KM = 6378.137
E_E    = 0.081819221456

# helper functions
## calculated denominator
def calc_denom(ecc, lat_rad):
  return math.sqrt(1.0-(ecc**2)*(math.sin(lat_rad)**2))

# initialize script arguments
o_x_km = float('nan') # ECEF x origin of SEZ frame (observer)
o_y_km = float('nan') # ECEF y origin of SEZ frame (observer)
o_z_km = float('nan') # ECEF z origin of SEZ frame (observer)
x_km = float('nan') # ECEF x position (object)
y_km = float('nan') # ECEF y position (object)
z_km = float('nan') # ECEF z position (object)

# parse script arguments
if len(sys.argv)==7:
    o_x_km = float(sys.argv[1])
    o_y_km = float(sys.argv[2])
    o_z_km  = float(sys.argv[3])
    x_km = float(sys.argv[4])
    y_km = float(sys.argv[5])
    z_km = float(sys.argv[6])
else:
    print(\
     'Usage: '\
     'python3 sez_to_ecef.py o_lat_deg o_lon_deg o_hae_km s_km e_km z_km'\
    )
    exit()


# write script below this line

# ecef to llh
# calculate longitude
lon_rad = math.atan2(o_y_km,o_x_km)
lon_deg = lon_rad*180.0/math.pi

# initialize lat_rad, r_lon_km, r_z_km
lat_rad = math.asin(o_z_km/math.sqrt(o_x_km**2+o_y_km**2+o_z_km**2))
r_lon_km = math.sqrt(o_x_km**2+o_y_km**2)
prev_lat_rad = float('nan')

# iteratively find latitude
c_E = float('nan')
count = 0
while (math.isnan(prev_lat_rad) or abs(lat_rad-prev_lat_rad)>10e-7) and count<5:
  denom = calc_denom(E_E,lat_rad)
  c_E = R_E_KM/denom
  prev_lat_rad = lat_rad
  lat_rad = math.atan((o_z_km+c_E*(E_E**2)*math.sin(lat_rad))/r_lon_km)
  count = count+1
lat_deg = lat_rad*180.0/math.pi 

lat_deg
# calculate hae
hae_km = r_lon_km/math.cos(lat_rad)-c_E

# ECEF Vector if station vs object is given
ecef_origin = [o_x_km, o_y_km, o_z_km]
ecef_position = [x_km, y_km, z_km]
r_sez = [[ecef_position[0]-ecef_origin[0]],
         [ecef_position[1]-ecef_origin[1]],
         [ecef_position[2]-ecef_origin[2]]]

# Inverse rotations
# rz_ry_inverse = (math.sin(lat_rad))
ry_inverse = [[math.sin(lat_rad),0,-math.cos(lat_rad)],
              [0,1,0],
              [math.cos(lat_rad), 0, math.sin(lat_rad)]]
rz_inverse = [[math.cos(lon_rad), math.sin(lon_rad), 0],
              [-math.sin(lon_rad), math.cos(lon_rad),0],
              [0,0,1]]

ry_rz_product = [
   [ry_inverse[0][0] * rz_inverse[0][0] + ry_inverse[0][1] * rz_inverse[1][0] + ry_inverse[0][2] * rz_inverse[2][0],
    ry_inverse[0][0] * rz_inverse[0][1] + ry_inverse[0][1] * rz_inverse[1][1] + ry_inverse[0][2] * rz_inverse[2][1],
    ry_inverse[0][0] * rz_inverse[0][2] + ry_inverse[0][1] * rz_inverse[1][2] + ry_inverse[0][2] * rz_inverse[2][2]],
    [ry_inverse[1][0] * rz_inverse[0][0] + ry_inverse[1][1] * rz_inverse[1][0] + ry_inverse[1][2] * rz_inverse[2][0],   
    ry_inverse[1][0] * rz_inverse[0][1] + ry_inverse[1][1] * rz_inverse[1][1] + ry_inverse[1][2] * rz_inverse[2][1],
    ry_inverse[1][0] * rz_inverse[0][2] + ry_inverse[1][1] * rz_inverse[1][2] + ry_inverse[1][2] * rz_inverse[2][2]],
    [ry_inverse[2][0] * rz_inverse[0][0] + ry_inverse[2][1] * rz_inverse[1][0] + ry_inverse[2][2] * rz_inverse[2][0],
    ry_inverse[2][0] * rz_inverse[0][1] + ry_inverse[2][1] * rz_inverse[1][1] + ry_inverse[2][2] * rz_inverse[2][1],
    ry_inverse[2][0] * rz_inverse[0][2] + ry_inverse[2][1] * rz_inverse[1][2] + ry_inverse[2][2] * rz_inverse[2][2]]
]

ry_rz_rsez_product = [
    [ry_rz_product[0][0] * r_sez[0][0] + ry_rz_product[0][1] * r_sez[1][0] + ry_rz_product[0][2] * r_sez[2][0]],
    [ry_rz_product[1][0] * r_sez[0][0] + ry_rz_product[1][1] * r_sez[1][0] + ry_rz_product[1][2] * r_sez[2][0]],
    [ry_rz_product[2][0] * r_sez[0][0] + ry_rz_product[2][1] * r_sez[1][0] + ry_rz_product[2][2] * r_sez[2][0]]
]

for row in ry_rz_rsez_product:
    print(row)
