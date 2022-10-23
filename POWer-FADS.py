from corrections_on_zw import *
from def_case import *
from int_z_minus_cst_dz import *
from make_graph import *
from Calc_V import *
from calc_k import *
from Txt_outputs import *
from interp_time_data import *
from Calc_zw_ini import *
from calc_k_rect_int___partial_unsat_soil_and_wet_well import *
import numpy as np
import random
import sys
import timeit

#read input folders as command line arguments
Input_folders =sys.argv[1:]
#Check if an input folder has been properly given
number_wells=len(Input_folders)
if(number_wells==0):
        print("ERROR")
        print ("POWeR-FADS.py needs as command line arg the name of the relevant input folder(s)")
        print("For example, run the following command : python POWeR-FADS.py Test_case_1")
        print("or : python POWeR-FADS.py Test_case_1 Test_case_2")
        sys.exit()
else:
        #Measure start time
        t0_start = timeit.default_timer()
        #State that the code begins
        if (number_wells==1):
                print("POWeR-FADS MODEL STARTS FOR 1 WELL: "+Input_folders[0])
        else:
                print("POWeR-FADS MODEL STARTS FOR "+str(number_wells)+" WELLS :")
                for counter,Input_folder in enumerate(Input_folders):
                        print(Input_folder)                        
        for counter,Input_folder in enumerate(Input_folders):
                print('***Simulation for '+Input_folder+" starts")
                #Read POWeR-FADS settings
                file_settings=np.loadtxt(fname=Input_folder+"/Inputs/settings.txt", dtype=float, usecols=1)
                write_fluxes=int(file_settings[0])
                display_graphical_outputs=int(file_settings[1])
                nbr_Interp_deltaT=int(file_settings[2])
                nbr_Interp_deltaT_MAX=int(file_settings[3])
                nbr_pts_rect_int_K  =int(file_settings[4])
                niter_max =int(file_settings[5])
                n_failmax =int(file_settings[6])
                epsilon =file_settings[7]
                Exchanges_well_matrix_through_well_bottom=file_settings[8]
                AveragingK=int(file_settings[9])
                calc_auto_zw0=int(file_settings[10])
                #If calc_auto_zw0=0 an initial value of zw0 must be defined
                if (calc_auto_zw0==0):
                        try:
                                value_zw0=float(file_settings[11])
                        except IndexError:
                                print("ERROR")
                                print("If calc_auto_zw0=0 you must add a user-defined value for zw0 on the next line in settings.txt")
                                print("for example, add this on the last line in settings.txt :")
                                print(" value_zw0(only_useful_if_calc_auto_zw0=0) 12 ")
                                sys.exit(1)
                else:
                        value_zw0=-99999.0
                #Read params_well
                ztube,zsurf_,zbot_f,zbot_m,rin_,rout_,rdrill_,wdrill_,deltaE_=np.loadtxt(fname=Input_folder+"/Inputs/params_well.txt", dtype=float, usecols=1)
                #Physical limits on zbot_f
                if (zbot_f<zbot_m):
                        zbot_f=zbot_m
                elif(zbot_f>zsurf_):
                        zbot_f=zsurf_
                #Read Time series of zm and zf
                file_zm_zf=params_well=np.loadtxt(fname=Input_folder+"/Inputs/zm_zf.txt", dtype=float,skiprows=1)
                times=file_zm_zf[:,0]
                vector_zf=file_zm_zf[:,2]
                vector_zm=file_zm_zf[:,1]
                nb_simulated_times=len(times)
                #Read params_mat and params_fract
                param_mat=np.loadtxt(fname=Input_folder+"/Inputs/params_matrix.txt", dtype=float, usecols=1)
                param_fract=np.loadtxt(fname=Input_folder+"/Inputs/params_fractures.txt", dtype=float, usecols=1)
                #interpolation of time series for zm and zf
                times_interp,nbr_dates_interp,periode_interp,vector_zf_interp,vector_zm_interp=interp_time_data(times,vector_zf,vector_zm,nbr_Interp_deltaT)
                #initialisation of zw 
                zw=Calc_zw_ini(calc_auto=calc_auto_zw0,zw0=value_zw0,zf_steady=vector_zf_interp[0],zm_steady=vector_zm_interp[0],zbot_f_=zbot_f,ztube_=ztube)
                #other initilaisations
                Zw_final,stockfract,stockmat,fluxmat,fluxfract=np.full((5,nbr_dates_interp),0.0)
                stockmat_foo,stockfract_foo,index_time=(0.0,0.0,-1)
                #START OF CALCULATIONS
                #cpu time counter start
                t1_start = timeit.default_timer()
                #time loop
                while (index_time<nbr_dates_interp-1):
                         index_time=index_time+1
                         zw_old_=zw
                         zf=0.5*vector_zf_interp[min(index_time+1,nbr_dates_interp-1)]+0.5*vector_zf_interp[index_time]
                         zm=0.5*vector_zm_interp[min(index_time+1,nbr_dates_interp-1)]+0.5*vector_zm_interp[index_time]
                         convergence=False
                         index_iterations=-1
                         zw_it=[-999999.0]*(niter_max+2)
                         Delta=[-999999.0]*(niter_max+2)
                         F=[-999999.0]*(niter_max+2)
                         J=[-999999.0]*(niter_max+2)
                         fail=0
                         #convergence loop
                         while (convergence==False):
                             index_iterations=index_iterations+1
                             #1st iteration : assume zw(n+1)=zw(n)
                             if (index_iterations==0):
                                     zw_it[index_iterations]=zw_old_
                             #following iterations : update the approximation of zw(n+1)
                             #NB : a random relaxaiton factor is used ot ease convergence
                             else:                 
                                 zw_it[index_iterations]=random.uniform(0.1, 1.0)*Delta[index_iterations-1]+zw_it[index_iterations-1] # random relaxation factor to help convegrence
                             #CALCULATION OF F
                             #find in what situation we are in for the matrix
                             case_m=def_case(zbot_m_=zbot_m,surface=zsurf_,zw__=zw_it[index_iterations],zi=zm)
                             #Calc Vm and its derivative
                             t4=calc_V(R_IN=rin_,R_OUT=rout_,R_DRILL=rdrill_,W_DRILL=wdrill_,flux_bottom=Exchanges_well_matrix_through_well_bottom,
                                case=case_m,parameters_medium=param_mat,_zi_=zm,z_bot_i=zbot_m,zsurf=zsurf_,deltaE=deltaE_,
                                nbr_points_interpolation_K=nbr_pts_rect_int_K ,zw_=zw_it[index_iterations],AveragingK_=AveragingK)
                             V_mat=t4[0]
                             dV_mat=t4[1]
                             #find in what situation we are in for the fractures
                             case_f=def_case(zbot_m_=zbot_f,surface=zsurf_,zw__=zw_it[index_iterations],zi=zf)
                             #Calc Vf and its derivative
                             t5=calc_V(R_IN=rin_,R_OUT=rout_,R_DRILL=rdrill_,W_DRILL=wdrill_,flux_bottom=False,case=case_f,parameters_medium=param_fract,_zi_=zf,z_bot_i=zbot_f,
                                   zsurf=zsurf_,deltaE=deltaE_,nbr_points_interpolation_K=nbr_pts_rect_int_K ,
                                   zw_=zw_it[index_iterations],AveragingK_=AveragingK)
                             Vfract=t5[0]
                             dVfract=t5[1]
                             #V is the sum of Vm and Vf
                             V_total=V_mat+Vfract
                             dV_total=dV_mat+dVfract
                             #conclude to calulate F
                             F[index_iterations]=zw_it[index_iterations]-zw_old_+periode_interp*V_total
                             #Is F a suitable approximation of 0 (i.e did we find a suitable approximation of zw(n+1) ?
                             #if yes :
                             if (abs(F[index_iterations])<epsilon):
                                 convergence=True
                                 #get value of zw
                                 zw_apriori=zw_it[index_iterations]
                                 zw=corrections_on_zw(zw_initial=zw_apriori,ztube_=ztube,zbot_m_=zbot_m)           
                                 Zw_final[index_time]=zw
                                 #Calc fluxes and sum of fluxes
                                 area=3.14159265359*rin_*rin_+wdrill_*3.14159265359*(rdrill_*rdrill_-rout_*rout_)
                                 fluxfract[index_time]=Vfract*area
                                 fluxmat[index_time]=V_mat*area
                                 stockfract_foo=stockfract_foo+fluxfract[index_time]*periode_interp
                                 stockfract[index_time]=stockfract_foo
                                 stockmat_foo=stockmat_foo+fluxmat[index_time]*periode_interp
                                 stockmat[index_time]=stockmat_foo   
                             #if no :
                             else:
                                #State that the simulation failed if too many iterations have already occured  
                                if (index_iterations>niter_max):
                                      index_iterations=-1
                                      zw_it=[-999999.0]*(niter_max+2)
                                      Delta=[-999999.0]*(niter_max+2)
                                      F=[-999999.0]*(niter_max+2)
                                      J=[-999999.0]*(niter_max+2)
                                      fail=fail+1
                                      #State the number of tries remaining
                                      if (fail==1):
                                              print("   Convergence failed at time "+"{:.2f}".format(times_interp[index_time])+" for 1 try, "+str(n_failmax-fail+1)+" remaining")
                                      else:
                                              print("   Convergence failed at time "+"{:.2f}".format(times_interp[index_time])+" for "+str(fail)+" trys, "+str(n_failmax-fail+1)+" remaining")
                                      #Possibility to start again from t=0 with a reduced time step if too many fails
                                      if (fail>n_failmax):
                                         print("  Convergence failed too many times")
                                         #Try to reduce time stem
                                         nbr_Interp_deltaT_old=nbr_Interp_deltaT
                                         nbr_Interp_deltaT=min(nbr_Interp_deltaT_old*2,nbr_Interp_deltaT_MAX)
                                         #If time step cannot be reduced, give up
                                         if (nbr_Interp_deltaT==nbr_Interp_deltaT_old):
                                            print("**No more reduction of time step is allowed : exiting, simulation failed")
                                            print
                                            sys.exit()
                                         #if we continue, w emust reinitialize everything with the new deltaT
                                         #interpolation of time series
                                         times_interp,nbr_dates_interp,periode_interp,vector_zf_interp,vector_zm_interp=interp_time_data(times,vector_zf,vector_zm,nbr_Interp_deltaT)
                                         print(str(f'  Retry with a reduced time step deltaT={periode_interp}'))
                                         #initialisation of zw
                                         zw=Calc_zw_ini(calc_auto=calc_auto_zw0,zw0=value_zw0,zf_steady=vector_zf_interp[0],zm_steady=vector_zm_interp[0],zbot_f_=zbot_f,ztube_=ztube)
                                         #other initilaisations
                                         Zw_final,stockfract,stockmat,fluxmat,fluxfract=np.full((5,nbr_dates_interp),0.0)
                                         stockmat_foo,stockfract_foo,index_time=(0.0,0.0,-1)
                                         break #This makes the whole simulation start again with a reduced time step
                                #Continue the convergence iterations if we did not reach the max number of iterations
                                else:
                                      #Find a new approximation of zw
                                      J[index_iterations]=1+periode_interp*dV_total
                                      Delta[index_iterations]=-F[index_iterations]/J[index_iterations]
                #cpu time counter end
                t1_stop = timeit.default_timer()
                #WRITE OUTPUTS
                Txt_outputs(path=Input_folder+'/Outputs',Zw_final_=Zw_final,stockmat_=stockmat,stockfract_=stockfract,
                    fluxmat_=fluxmat,fluxfract_=fluxfract,times_interp_=times_interp,write_fluxes_=write_fluxes)                      
                if (display_graphical_outputs!=0):
                    make_graph(outputfile=Input_folder+"/Outputs/Output_graph.png",
                           zm=vector_zm_interp,zf=vector_zf_interp,zw=Zw_final,
                           time=times_interp,flux_m=3600*fluxmat,flux_f=-3600*fluxfract,zbot_f_=zbot_f)
                print("   Simulation Successful for "+Input_folder+" in "+"{:.2f}".format(t1_stop-t1_start)+" seconds of CPU time")
t0_stop = timeit.default_timer()
print("ALL WELLS WERE SIMULATED SUCCESSFULLY")
print("Total ellapsed time : "+"{:.2f}".format(t0_stop-t0_start)+" seconds")
