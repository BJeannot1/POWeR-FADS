from calc_k import *
from calc_k_rect_int___partial_unsat_soil_and_wet_well import *
from int_z_minus_cst_dz import *
#Calculate Vf or Vm
def calc_V(flux_bottom,case,parameters_medium,_zi_,z_bot_i,zw_,zsurf,R_IN,R_OUT,R_DRILL,W_DRILL,deltaE,nbr_points_interpolation_K,AveragingK_):
    Ks_loc=parameters_medium[0]
    V_medium=0.0
    dV_medium=0.0
    if (case==1234 or case==2134 or case==1243 or case==2314 or case==2143):
         t1=calc_k_rect_int___partial_unsat_soil_and_wet_well(
             zbot_i=z_bot_i,
             ztop=max(z_bot_i,min(zw_,zsurf)),
             zi_=_zi_,
             zw_piezo=zw_,
             ndeltz=nbr_points_interpolation_K,
             parameters_mil=parameters_medium,
             avgK=AveragingK_)      
         K_INT=t1[0]
         dK_INT=t1[1]         
         vint=(2*R_DRILL*(zw_-_zi_)*K_INT)/((R_IN*R_IN+W_DRILL*R_DRILL*R_DRILL-W_DRILL*R_OUT*R_OUT)*deltaE)
         u=(zw_-_zi_)
         uprime=1.0
         v=K_INT
         vprime=dK_INT
         dvint=(2*R_DRILL/((R_IN*R_IN+W_DRILL*R_DRILL*R_DRILL-W_DRILL*R_OUT*R_OUT)*deltaE))*((u*vprime)+(v*uprime))
         if (flux_bottom==True):
             vflux_bottom=Ks_loc*(zw_-_zi_)/deltaE #seulement pour la matrice
             dvflux_bottom=Ks_loc/deltaE
         else:
             vflux_bottom=0
             dvflux_bottom=0
         V_medium=V_medium+vflux_bottom+vint
         dV_medium=dV_medium+dvflux_bottom+dvint
         if (V_medium<0 or vflux_bottom*vint<0):
             print("ERROR")            
             exit()
    elif (case==1324 or case==3124 or case==1342 or case==3214 or case==3142):
         t2=int_z_minus_cst_dz(zw_piezo=zw_,limit_inf=min(zsurf,max(zw_,z_bot_i)),limit_sup=max(z_bot_i,min(_zi_,zsurf)),cst=_zi_)
         integrale=t2[0]
         dintegrale=t2[1]
         v_piezo_vide=(2*R_DRILL*Ks_loc*integrale)/((R_IN*R_IN+W_DRILL*R_DRILL*R_DRILL-W_DRILL*R_OUT*R_OUT)*deltaE)
         dv_piezo_vide=(2*R_DRILL*Ks_loc*dintegrale)/((R_IN*R_IN+W_DRILL*R_DRILL*R_DRILL-W_DRILL*R_OUT*R_OUT)*deltaE)
         t3=calc_k_rect_int___partial_unsat_soil_and_wet_well(
             zbot_i=z_bot_i,
             ztop=max(z_bot_i,min(zw_,zsurf)),
             zi_=_zi_,
             zw_piezo=zw_,
             ndeltz=nbr_points_interpolation_K,
             parameters_mil=parameters_medium,
             avgK=AveragingK_)
         K_INT=t3[0]
         dK_INT=t3[1] 
         v_piezo_plein=2*R_DRILL*(zw_-_zi_)*K_INT/((R_IN*R_IN+W_DRILL*R_DRILL*R_DRILL-W_DRILL*R_OUT*R_OUT)*deltaE)
         u=(zw_-_zi_)
         uprime=1.0
         v=K_INT
         vprime=dK_INT
         dv_piezo_plein=(2*R_DRILL/((R_IN*R_IN+W_DRILL*R_DRILL*R_DRILL-W_DRILL*R_OUT*R_OUT)*deltaE))*((u*vprime)+(v*uprime))
         if (flux_bottom==True):
             v_flux_bottom=Ks_loc*(zw_-_zi_)/deltaE 
             dv_flux_bottom=Ks_loc/deltaE
         else:
             v_flux_bottom=0.0
             dv_flux_bottom=0.0
         V_medium=V_medium+v_piezo_vide+v_piezo_plein+v_flux_bottom
         dV_medium=dV_medium+dv_piezo_vide+dv_piezo_plein+dv_flux_bottom
         if (V_medium>0 or v_piezo_vide*v_piezo_plein<0 or v_piezo_plein*v_flux_bottom<0 ):          
             print("ERROR")            
             exit()
    elif (case==1432 or case==1423):
         V_medium=0.0
         dV_medium=0.0
    return V_medium,dV_medium
