#Function used ot initialize zw. based on te hypothesis of steady state zf and zm
def Calc_zw_ini(calc_auto,zw0,zf_steady,zm_steady,zbot_f_,ztube_):
    if (calc_auto!=0):
        if (zf_steady>=zbot_f_):
            zw=min(zf_steady,ztube_)
        else:
            zw = min(zbot_f_,zm_steady)
    else:
        zw = zw0
    return zw
