#fFunction to find what is the order of zsurf,zi,zw,zbot_m
def def_case(zbot_m_,surface,zw__,zi):
    if ((surface>=zw__)and (zw__>=zi)and (zi>=zbot_m_)):
        case=1234
    elif((surface>=zi)and (zi>=zw__)and (zw__>=zbot_m_)):
        case=1324
    elif((zw__>=surface)and (surface>=zi)and (zi>=zbot_m_)):
        case=2134
    elif((zi>=surface)and (surface>=zw__)and (zw__>=zbot_m_)):
        case=3124
    elif((surface>=zi)and (zi>=zbot_m_)and (zbot_m_>=zw__)):
        case=1342
    elif((surface>=zw__)and (zw__>=zbot_m_)and (zbot_m_>=zi)):
        case=1243        
    elif((surface>=zbot_m_)and (zbot_m_>=zi)and (zi>=zw__)):
        case=1432
    elif((surface>=zbot_m_)and (zbot_m_>=zw__)and (zw__>=zi)):
        case=1423
    elif((zw__>=zi)and (zi>=surface)and (surface>=zbot_m_)):
        case=2314
    elif((zi>=zw__)and (zw__>=surface)and (surface>=zbot_m_)):
        case=3214
    elif((zw__>=surface)and (surface>=zbot_m_)and (zbot_m_>=zi)):
        case=2143
    elif((zi>=surface)and (surface>=zbot_m_)and (zbot_m_>=zw__)):
        case=3142
    else:
        print ("ERROR")       
        exit()
    return case
