#calculate the int of (z-cste) for z varying between limit_inf and limit_sup
#Also calculate its derivative
def int_z_minus_cst_dz(limit_inf,limit_sup,cst,zw_piezo):
    if (limit_inf==limit_sup):
        Int=0.0
        dInt=0.0
    elif(limit_sup==cst):
        Int=-0.5*(cst-limit_inf)**2.0
        if (limit_inf!=zw_piezo):
            dInt=0.0
        else:
            dInt=cst-zw_piezo  
    else:
        Int=0.5*limit_sup*limit_sup-cst*limit_sup-0.5*limit_inf*limit_inf+limit_inf*cst
        if (limit_inf!=zw_piezo):
            dInt=0.0
        else:
            dInt=cst-zw_piezo            
    return Int,dInt
