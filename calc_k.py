#Mualem Van Genuchten to calculate unsaturated K and its derivative
def calc_k(pressure,parameters):
    #
    Ks=parameters[0]
    ngenuch=parameters[1]    
    alpha=parameters[2]
    k_min=parameters[3]
    
    if (pressure>-0.001):
            K=Ks
            dK=0.0
    else:
                          l_=0.5
                          t1 = pressure
                          t2 = abs(t1)
                          t4 = (alpha * t2) ** ngenuch
                          t5 = 1.0 + t4
                          t7 = 1.0 - 1.0 / ngenuch
                          t8 = t5 ** t7
                          SE = 1.0 / t8
                          t10 = SE ** l_
                          t12 = ngenuch - 1.0
                          t15 = SE ** (ngenuch / t12)
                          t16 = 1.0- t15
                          t17 = t16 ** t7
                          t18 = 1.0 - t17
                          t19 = t18 ** 2.0
                          K = max(k_min,Ks * t10 * t19)
                          t20=abs(t1)/t1
                          t26=-t16                          
                          dK=t12*t20*(((l_+2)*t15-l_)*t17-l_*t26)*t18*Ks*t10*t4/t2/t5/t26
    return K,dK
