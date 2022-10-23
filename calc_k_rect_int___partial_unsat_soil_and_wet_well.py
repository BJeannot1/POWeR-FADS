#integration of K in ndelt parts from zbot_i to ztop, in a context where
#the well is full of water between zbot_i and ztop
from calc_k import *
def calc_k_rect_int___partial_unsat_soil_and_wet_well(zbot_i,ztop,zi_,zw_piezo,ndeltz,parameters_mil,avgK):
    if (zbot_i==ztop):
        K_int=0.0
        dK_int=0.0
    else:
        Kslocal=parameters_mil[0]
        pressure_sol_top = zi_ - ztop
        K_int=0.0
        dK_int=0.0
        #step 1 : check if this subroutine has been called innide its domain of application
        if (zw_piezo<ztop):
            print("error : the well is not full of water between zbot_i and ztop")
            exit()
        #Case 1 : Totally unsaturated soil
        if (zi_ <= zbot_i):
            t100=ztop-zbot_i
            deltz = t100 / ndeltz
            prdelt = pressure_sol_top + (0.5 * deltz+t100)
            for j21 in range(1, ndeltz+1,1):
                prdelt = prdelt - deltz
                t61=calc_k(pressure=prdelt,parameters=parameters_mil)
                ksol=t61[0]
                dksol=t61[1]
                kpiezo=Kslocal
                if (avgK==0.0): #arithmetic mean
                    moy_k=(ksol+kpiezo)*0.5 
                    dmoy_k=(dksol)*0.5 
                elif (avgK==1.0): #geometric mean 
                    moy_k=(ksol*kpiezo)**0.5 
                    dmoy_k=kpiezo*dksol*0.5*(ksol*kpiezo)**(-0.5) 
                else:
                    print("ERROR")
                    exit()
                K_int=K_int+ moy_k*deltz
                dK_int=dK_int+ ((0.5-j21)*dmoy_k/ndeltz)*deltz+moy_k/ndeltz
        #Case 2 : Partly unsaturated soil
        elif (zi_ <= ztop):
            t100=ztop-zi_
            deltz = t100 / ndeltz
            prdelt = pressure_sol_top + (0.5 * deltz+t100)
            for j21 in range(1, ndeltz+1,1):
                prdelt = prdelt - deltz
                t61=calc_k(pressure=prdelt,parameters=parameters_mil)
                ksol=t61[0]
                dksol=t61[1]
                kpiezo=Kslocal
                if (avgK==0.0): #arithmetic mean
                    moy_k=(ksol+kpiezo)*0.5 
                    dmoy_k=(dksol)*0.5 
                elif (avgK==1.0): #geometric mean 
                    moy_k=(ksol*kpiezo)**0.5 
                    dmoy_k=kpiezo*dksol*0.5*(ksol*kpiezo)**(-0.5) 
                else:
                    print("ERROR")
                    exit()
                K_int=K_int+ moy_k*deltz
                dK_int=dK_int+ ((0.5-j21)*dmoy_k/ndeltz)*deltz+moy_k/ndeltz
            K_int=K_int+(zi_-zbot_i)*Kslocal
        #Case 3 : Saturated soil
        else:
            K_int=K_int+(ztop-zbot_i)*Kslocal            
            if (ztop!=zw_piezo):
                dK_int=0.0
            else:
                dK_int=Kslocal        
    return K_int,dK_int
