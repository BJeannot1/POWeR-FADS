import numpy as np
#Function which inteprolates time, zf and zm linearly
def interp_time_data(temps_,vector_zf_,vector_zm_,nbr_Interp_deltaT_):
    nbr_dates=len(temps_)
    temps_interp=np.linspace(temps_[0],temps_[-1],num=(nbr_dates-1)*nbr_Interp_deltaT_+1)
    vector_zf_interp=np.interp(x=temps_interp,xp=temps_,fp=vector_zf_)
    vector_zm_interp=np.interp(x=temps_interp,xp=temps_,fp=vector_zm_)
    nbr_dates_interp=len(temps_interp)
    periode_interp=(temps_[-1]-temps_[0])/(nbr_dates_interp-1)
    return temps_interp,nbr_dates_interp,periode_interp,vector_zf_interp,vector_zm_interp
