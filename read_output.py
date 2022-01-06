import numpy as np 
import struct
import matplotlib.pyplot as plt 
from matplotlib import cm
import sys
import glob
from matplotlib.gridspec import GridSpec
from matplotlib.colors import TwoSlopeNorm
import array
import os
from scipy.constants import k,mu_0
from matplotlib.collections import LineCollection

RS = 6.96342e8 # meter
MP = 1.67262192369e-27 # kg, mass of proton
ME = 9.10938356E-31 # kg, mass of electron
KB = 1.380649e-23 # Boltzmann constant

def read_grid():
    # read nx, xgrid-------------------------
    file_grid = open('./output/xgrid.dat', 'rb')

    nx = int((struct.unpack("d",file_grid.read(8)))[0])
    xgrid = np.zeros(nx)

    for i in range(nx):
        xgrid[i] = (struct.unpack("d",file_grid.read(8)))[0]

    file_grid.close()

    return xgrid

def read_1d_arr(filename,offset=8):
    file_ = open(filename,'rb')

    nx = int((struct.unpack("d",file_.read(offset)))[0])

    data = np.fromfile(filename,offset=offset)
    return nx,data

def read_output(filename, nvar,nx):
    arr = array.array('d')
    arr.fromfile(open(filename, 'rb'), nvar*nx + 1)

    t = arr[0]
    uu = np.copy(np.reshape(np.asarray(arr[1:len(arr)]), [nx,nvar], order='C'))
    arr = -1
    return t, uu

def read_parameter_from_input(kw):
    file = open('input','r')

    info = file.readlines()

    for i in range(len(info)):
        if info[i].find(kw)!=-1:
            str_to_read = info[i]

            value = float(str_to_read[0:(str_to_read.find(';')-1)])

            return value 
    
    print('The keyword is not found !!!')
    return None

offset = 8 

nvar = 6
nvarPr = 3

adiabaticIdx_I = read_parameter_from_input('adiabatic_index Ion')
print('Adiabatic Index Ion = {:.5f}'.format(adiabaticIdx_I))

adiabaticIdx_E = read_parameter_from_input('adiabatic_index Electron')
print('Adiabatic Index Ion = {:.5f}'.format(adiabaticIdx_E))

xgrid = read_grid()
nx = len(xgrid)

print('nx = ', nx)
print('x0 = {:.3f}, x1 = {:.3f}'.format(xgrid[0]/RS,
    xgrid[-1]/RS))


xgrid_norm = xgrid/RS 


nx,Br = read_1d_arr('./output/Br.dat',offset)
nx,fExp = read_1d_arr('./output/fExp.dat',offset)
# plt.plot(xgrid_norm, Br)
# plt.plot(xcellinterface/RS,BrCI,'C1--')
# plt.show()

print('Br(1AU) = {:.3f} nT'.format(Br[-1]/1e-9))


x_to_cut = 215
ind_cut = np.where(np.abs(xgrid_norm - x_to_cut) == 
    np.min(np.abs(xgrid_norm - x_to_cut)))[0][0]


files = sorted(glob.glob('output/out*dat'))
nfiles = len(files)

xlim = [-5,xgrid_norm.max()]

segs_u = np.zeros([nfiles, nx, 2])
segs_n = np.zeros([nfiles, nx, 2])
segs_Ti = np.zeros([nfiles, nx, 2])
segs_Te = np.zeros([nfiles, nx, 2])

time = np.zeros(nfiles)

u_cut = np.zeros(nfiles)
n_cut = np.zeros(nfiles)
Ti_cut = np.zeros(nfiles)
Te_cut = np.zeros(nfiles)

for nt in range(nfiles):
# nt = nfiles-1
    filename = files[nt]
    t, uu = read_output(filename, nvar, nx)

    print('Read file: t = {:.3f}'.format(t))


    rho = uu[:,0]
    u_ = uu[:,1]
    pi_ = uu[:,2]
    pe_ = uu[:,3]
    eOut = uu[:,4]
    eIn = uu[:,5]

    n = rho/(MP+ME)

    zout_ = np.sqrt(4*eOut/rho)
    zin_ = np.sqrt(4*eIn/rho)

    Ti =  pi_ / n / k 
    Te = pe_ / n / k

    Va = Br/np.sqrt(rho * mu_0)
    Cs = np.sqrt((adiabaticIdx_I * pi_ + adiabaticIdx_E * pe_)/rho)

    # print(max(Va)/1e3, max(u_)/1e3, max(Cs)/1e3)
    # print(min(Va)/1e3, min(u_)/1e3, min(Cs)/1e3)
    # continue

    print(n[0]/1e16,Ti[0]/1e5,Te[0]/1e5,u_[0]/1e3,zout_[0]/1e3,n[-1]/1e6)

    u_cut[nt] = u_[ind_cut]
    n_cut[nt] = n[ind_cut]
    Ti_cut[nt] = Ti[ind_cut]
    Te_cut[nt] = Te[ind_cut]

    segs_n[nt, :, 1] = n[:] 
    segs_u[nt, :, 1] = u_[:]
    segs_Ti[nt, :, 1] = Ti[:] 
    segs_Te[nt, :, 1] = Te[:]

    segs_n[nt, :, 0] = xgrid_norm[:]
    segs_u[nt, :, 0] = xgrid_norm[:]
    segs_Ti[nt, :, 0] = xgrid_norm[:]
    segs_Te[nt, :, 0] = xgrid_norm[:]
    time[nt] = t

    # if nt==nfiles-1:
    #     plt.plot(xgrid_norm, u_, marker='.', ls='none')
    #     plt.show()

    #     sys.exit()

    if len(glob.glob('./figure/{:04d}.png'.format(nt)))>0:
        continue

    fig = plt.figure(figsize=[7,12])

    sub = fig.add_subplot(311)
    sub.plot(xgrid_norm,n,label=r'$N$/$m^{-3}$')
    sub.set_ylabel(r'$N$/$m^{-3}$')
    sub.set_xlim([1.0,215])
    sub.set_xscale('log', base = 10)
    sub.set_yscale('log', base = 10)

    sub_t = sub.twinx()
    sub_t.plot(xgrid_norm,Ti,color='C1',label=r'$T_i$')
    sub_t.plot(xgrid_norm,Te,color='C2',label=r'$T_e$')
    sub_t.legend()
    sub.set_xscale('log', base = 10)
    sub.set_yscale('log', base = 10)

    sub = fig.add_subplot(312)
    sub.plot(xgrid_norm,u_,label=r'$U$')
    sub.plot(xgrid_norm,Va,label=r'$V_{A}$')
    sub.plot(xgrid_norm,Cs,label=r'$C_{s}$')
    sub.legend()
    sub.set_xlim(xlim)
    sub.set_ylim([-100,800e3])

    sub = fig.add_subplot(313)
    sub.plot(xgrid_norm,zout_,label=r'$Z_{o}$')
    sub.plot(xgrid_norm,zin_,label=r'$Z_{i}$')
    sub.legend()
    sub.set_xlim(xlim)

    fig.suptitle(r'$t={:.3f}$'.format(t))

    # plt.show()
    fig.savefig('./figure/{:04d}.png'.format(nt))
    plt.close(fig)




fig = plt.figure(figsize=[10,14])

sub = fig.add_subplot(411)
lc_n = LineCollection(segs_n, array=time, cmap=plt.cm.plasma)
sub.add_collection(lc_n)
sub.set_ylabel(r'$n$ m$^{-3}$', fontsize=14)
sub.set_xlim([1,215])
# sub.set_ylim([np.min(segs_n[:,:,1]), np.max(segs_n[:,:,1])])
sub.set_label(r'$x$')
sub.set_xscale('log',base=10)
sub.set_yscale('log',base=10)

sub0 = sub

sub = fig.add_subplot(412)
lc_u = LineCollection(segs_u, array=time, cmap=plt.cm.plasma)
sub.add_collection(lc_u)
sub.set_ylabel(r'$U$ m/s', fontsize=14)
sub.set_xlim(xlim)
# sub.set_ylim([np.min(segs_ux[:,:,1]), np.max(segs_ux[:,:,1])])
sub.set_ylim([-1e3,500e3])

sub1 = sub

sub = fig.add_subplot(413)
lc_Ti = LineCollection(segs_Ti, array=time, cmap=plt.cm.plasma)
sub.add_collection(lc_Ti)
sub.set_ylabel(r'$T_i$ K', fontsize=14)
sub.set_xlim([1,215])
sub.set_label(r'$x$')
sub.set_xscale('log',base=10)
sub.set_yscale('log',base=10)

sub2 = sub

sub = fig.add_subplot(414)
lc_Te = LineCollection(segs_Te, array=time, cmap=plt.cm.plasma)
sub.add_collection(lc_Te)
sub.set_ylabel(r'$T_e$ K', fontsize=14)
sub.set_xlim([1,215])
sub.set_label(r'$x$')
sub.set_xscale('log',base=10)
sub.set_yscale('log',base=10)

cb = fig.colorbar(lc_u, ax=[sub0, sub1, sub2, sub])
cb.ax.set_title('time/s',fontsize=12)

fig.savefig('./plot_collective.png', dpi=400)
plt.close(fig)



fig = plt.figure(figsize=[8,10])
sub = fig.add_subplot(311)
sub.plot(time,n_cut)
sub.set_ylabel(r'$n$ m$^{-3}$')
sub.tick_params(axis='x', which='both', direction='in')

sub = fig.add_subplot(312)
sub.plot(time,u_cut)
sub.set_ylabel(r'$V$ m/s')
sub.tick_params(axis='x', which='both', direction='in')

sub = fig.add_subplot(313)
sub.plot(time, Ti_cut)
sub.set_ylabel(r'$T_i$ K')
sub.set_xlabel(r'$t$ s')
sub.tick_params(axis='y', which='both', colors='C0')

subt = sub.twinx()
subt.plot(time, Te_cut, color='C1')
subt.set_ylabel(r'$T_e$ K')
subt.tick_params(axis='y', which='both', colors='C1')

fig.suptitle('$r = {:.3f}$'.format(x_to_cut), fontsize=14)
fig.subplots_adjust(hspace=0,top=0.95)

fig.savefig('./time_evolution.png',dpi=400)
plt.close(fig)
