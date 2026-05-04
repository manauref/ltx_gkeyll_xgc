import numpy as np
import h5py 
import adios2 as ad
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
from matplotlib.colors import LogNorm

#datadir = "/home/george/proj/ltx/turb/lowres-d2/"
datadir = "/home/george/proj/ltx/"
step = 9000
step_0 = 100

f2dfile = datadir+"xgc.f2d.%05d.bp"%(step)

with ad.Stream(datadir+"xgc.mesh.bp","rra") as f:
    rz = f.read("rz")
    conn = f.read("nd_connect_list")
    psi_surf = f.read("psi_surf")
    surf_idx = f.read("surf_idx") -1
    surf_len = f.read("surf_len")

with ad.Stream(datadir+"xgc.equil.bp","rra") as f:
    x_psi = f.read("eq_x_psi")
    x_r = f.read("eq_x_r")
    x_z = f.read("eq_x_z")
    axis_r = f.read("eq_axis_r")
    axis_z = f.read("eq_axis_z")

with ad.Stream(f2dfile,"rra") as f:
#    ne = np.average(f.read("e_den"),axis=0)
#    ne = f.read("e_den_en")[:,0]
    ne = f.read("e_den")[:]
    ni = f.read("i_den")[:]
    ui = f.read("i_u_para")[:]*1.0e-3
    ue = f.read("e_u_para")[:]

#    Teperp = np.average(f.read("e_T_perp"),axis=0)
#    Tepara = np.average(f.read("e_T_para"),axis=0)
#    Tiperp = np.average(f.read("i_T_perp"),axis=0)
#    Tipara = np.average(f.read("i_T_para"),axis=0)
    Teperp = f.read("e_T_perp")
    Tepara = f.read("e_T_para")
    Tiperp = f.read("i_T_perp")
    Tipara = f.read("i_T_para")
    Te = (2.0*Teperp + Tepara)/3.0
    Ti = (2.0*Tiperp + Tipara)/3.0

f2dfile_init = datadir+"xgc.f2d.%05d.bp"%(step_0)
with ad.Stream(f2dfile_init,"rra") as f:
#    ne_0 = np.average(f.read("e_den"),axis=0)
#    ne_0 = f.read("e_den_en")[:,0]
    ne_0 = f.read("e_den")[:]
    ui_0 = f.read("i_u_para")[:]*1.0e-3
    Teperp_0 = f.read("e_T_perp")
    Tepara_0 = f.read("e_T_para")
    Tiperp_0 = f.read("i_T_perp")
    Tipara_0 = f.read("i_T_para")
    Te_0 = (2.0*Teperp_0 + Tepara_0)/3.0
    Ti_0 = (2.0*Tiperp_0 + Tipara_0)/3.0

nsurf = len(psi_surf)-1
ne_omp = np.zeros(nsurf)
Te_omp = np.zeros(nsurf)
Ti_omp = np.zeros(nsurf)
ne_omp_0 = np.zeros(nsurf)
Te_omp_0 = np.zeros(nsurf)
Ti_omp_0 = np.zeros(nsurf)
psiN = psi_surf[1:]/x_psi
for i in range(1,nsurf+1):
    surf_idx_local = surf_idx[i,0:surf_len[i]]
#    idx_outer = surf_idx[i,rz[ surf_idx[i,0:surf_len[i]],0] > axis_r]

    omp_idx = surf_idx[i,np.argmax(rz[surf_idx_local,0] )]
    ne_omp[i-1] = ne[omp_idx]
    Te_omp[i-1] = Te[omp_idx]
    Ti_omp[i-1] = Ti[omp_idx]
    ne_omp_0[i-1] = ne_0[omp_idx]
    Te_omp_0[i-1] = Te_0[omp_idx]
    Ti_omp_0[i-1] = Ti_0[omp_idx]

#plt.figure(figsize=[5,8])
plt.subplot(311)
plt.plot(psiN,ne_omp)
plt.plot(psiN,ne_omp_0)
plt.subplot(312)
plt.plot(psiN,Te_omp)
plt.plot(psiN,Te_omp_0)
plt.subplot(313)
plt.plot(psiN,Ti_omp)
plt.plot(psiN,Ti_omp_0)
plt.savefig("moments.png")

def write_grid(filename,rz,conn):
    f = h5py.File(filename,"w")
    f.create_dataset("subplot00_nodes_xvalues",data=rz[:,0])
    f.create_dataset("subplot00_nodes_yvalues",data=rz[:,1])
    f.create_dataset("subplot00_line_connections",data=conn)
    f.close()

def write_moms(filename,x,y,z):
    f = h5py.File(filename,"w")
    f.create_dataset("subplot00_xvalues",data=x)
    f.create_dataset("subplot00_yvalues",data=y)
    f.create_dataset("subplot00_zvalues",data=z)
    f.close()



def write_omp(filename,psi,ne,Te,Ti):
    f = h5py.File(filename,"w")
    f.create_dataset("subplot00_line0_xvalues",data=psi)
    f.create_dataset("subplot00_line0_yvalues",data=ne)
    f.create_dataset("subplot01_line0_xvalues",data=psi)
    f.create_dataset("subplot01_line0_yvalues",data=Te)
    f.create_dataset("subplot02_line0_xvalues",data=psi)
    f.create_dataset("subplot02_line0_yvalues",data=Ti)
    f.close()

h5filename_final = "ltx103955-04_v0_xgc_final_moments_omp.h5"
h5filename_init = "ltx103955-04_v0_xgc_init_moments_omp.h5"
write_omp(h5filename_final,psiN,ne_omp,Te_omp,Ti_omp)
write_omp(h5filename_init,psiN,ne_omp_0,Te_omp_0,Ti_omp_0)

rz = ad.Stream(datadir+"xgc.mesh.bp","rra").read("rz")
conn = ad.Stream(datadir+"xgc.mesh.bp","rra").read("nd_connect_list")
write_grid("ltx_xgc_li863mg_grid_RZ.h5",rz,conn)

write_moms("ltx_xgc_li863mg_final_elc_den_RZ.h5",rz[:,0],rz[:,1],ne)
write_moms("ltx_xgc_li863mg_final_elc_upar_RZ.h5",rz[:,0],rz[:,1],ue)
write_moms("ltx_xgc_li863mg_final_elc_tpar_RZ.h5",rz[:,0],rz[:,1],Tepara)
write_moms("ltx_xgc_li863mg_final_elc_tperp_RZ.h5",rz[:,0],rz[:,1],Teperp)

write_moms("ltx_xgc_li863mg_final_ion_den_RZ.h5",rz[:,0],rz[:,1],ni)
write_moms("ltx_xgc_li863mg_final_ion_upar_RZ.h5",rz[:,0],rz[:,1],ui)
write_moms("ltx_xgc_li863mg_final_ion_tpar_RZ.h5",rz[:,0],rz[:,1],Tipara)
write_moms("ltx_xgc_li863mg_final_ion_tperp_RZ.h5",rz[:,0],rz[:,1],Tiperp)
                    


