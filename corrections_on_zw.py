def corrections_on_zw(zw_initial,ztube_,zbot_m_):
     #limits on Zw
     # ztube>=zw>=zbot_m_m  
     if (zw_initial>ztube_):
         zw_corrected=ztube_
     elif (zw_initial<zbot_m_):
         zw_corrected=zbot_m_
     else:
         zw_corrected=zw_initial
     return zw_corrected
