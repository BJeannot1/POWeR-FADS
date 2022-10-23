import os
from emptyfolder import *
def Txt_outputs(path,Zw_final_,stockmat_,stockfract_,fluxmat_,fluxfract_,times_interp_,write_fluxes_):
        #Check if output directory exists
        isExist = os.path.exists(path)
        #if not, create it
        if (isExist==False):
                os.mkdir(path)
        #if yes, empty it
        else:
                emptyfolder(path)
        #Open output files
        h10=open(path+"/zw.txt","w")
        h10.write("time(s),zw(m)\n")
        if (write_fluxes_==1):
                h20=open(path+"/Sum_fluxes.txt","w")
                h20.write("time(s),Sum_fluxes_from_well_towards_mat(m3),Sum_fluxes_from_well_towards_fract(m3)\n")
                h21=open(path+"/Fluxes.txt","w")
                h21.write("time(s),flux_from_well_towards_mat(m3/s),flux_from_well_towards_fract(m3/s)\n")
        for t in range(0,len(times_interp_),1):
                 h10.write(str(f' {times_interp_[t]},{Zw_final_[t]:.11f}\n'))
                 if (write_fluxes_==1):
                        h20.write(str(f' {times_interp_[t]},{stockmat_[t]:.11f},{stockfract_[t]:.11f}\n'))
                        h21.write(str(f' {times_interp_[t]},{fluxmat_[t]:.11f},{fluxfract_[t]:.11f}\n'))
        h10.write("\n")
        h10.close()
        if (write_fluxes_==1):
                h20.write("\n")
                h21.write("\n")     
                h20.close()
                h21.close()
