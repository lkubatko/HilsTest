import numpy as np
import random
import math
#import matplotlib.pyplot as plt

########################################################################
## Simulate gene trees from a species tree under the coalescent model ##
########################################################################


##### Number of trees to simulate
ntrees = 50000

##### Set species trees parameters

### Option 2: asymmetric species tree

thetaA = 0.005
thetaB = 0.01
thetaC = 0.005
thetaD = 0.005
thetaAB = 0.005
thetaABC = 0.005
thetaABCD = 0.005

# taus in coalescent units = number of generations per 2N
tau1 = 0.25
tau2 = 0.25
tau3 = 0.5

count = np.zeros(34)

for i in range(ntrees):
    tAB = np.random.exponential(1.0)
    if tAB<tau2:   # A and B coalesce in their ancestral interval
        tApr = str(tau1*thetaA+tAB*thetaAB)
        tBpr = str(tau1*thetaB+tAB*thetaAB)
        tABC = np.random.exponential(1.0)
        if tABC<tau3:  # AB and C coalesce in their ancestral interval
            tABpr = str(thetaAB*(tau2-tAB)+thetaABC*tABC)
            tABCD = np.random.exponential(1.0) # time for ABC to coalesce with D
            tABCpr = str(thetaABC*(tau3-tABC)+thetaABCD*tABCD)
            tCpr = str(thetaC*(tau1+tau2)+thetaABC*tABC)
            tDpr = str(thetaD*(tau1+tau2+tau3)+thetaABCD*tABCD)
            print '[1](((3:'+tApr+',4:'+tBpr+'):'+tABpr+',2:'+tCpr+'):'+tABCpr+',1:'+tDpr+');'
            count[0] = count[0]+1
            
        else: # AB and C do not coalesce in their ancestral interval
            ran1 = np.random.uniform()
            if ran1<1.0/3.0:  #AB and C coalesce first, above the root
                tABC = np.random.exponential(1.0/3.0)
                tABCD = np.random.exponential(1.0)
                tABpr = str(thetaAB*(tau2-tAB)+thetaABC*tau3+thetaABCD*tABC)
                tABCpr = str(thetaABCD*tABCD)
                tCpr = str(thetaC*(tau1+tau2+tau3)+thetaABCD*tABC)
                tDpr = str(thetaD*(tau1+tau2+tau3)+thetaABCD*tABCD)
                print '[1](((3:'+tApr+',4:'+tBpr+'):'+tABpr+',2:'+tCpr+'):'+tABCpr+',1:'+tDpr+');'
                count[1] = count[1]+1

            elif ran1<2.0/3.0: # AB and D coalesce first, above the root
                tABD = np.random.exponential(1.0/3.0)
                tABpr = str(thetaAB*(tau2-tAB)+thetaABC*tau3+thetaABCD*tABD)
                tDpr = str(thetaD*(tau1+tau2+tau3)+thetaABCD*tABD)
                tABCD = np.random.exponential(1.0)
                tABDpr = str(thetaABCD*tABCD)
                tCpr = str(thetaC*(tau1+tau2)+thetaABC*tau3+thetaABCD*(tABD+tABCD))
                print '[1](((3:'+tApr+',4:'+tBpr+'):'+tABpr+',1:'+tDpr+'):'+tABDpr+',2:'+tCpr+');'
                count[26] = count[26] + 1
                                            
            else: # C and D coalesce first, above the root
                tCD = np.random.exponential(1.0/3.0)
                tCpr = str(thetaC*(tau1+tau2)+thetaABC*tau3+thetaABCD*tCD)
                tDpr = str(thetaD*(tau1+tau2+tau3)+thetaABCD*tCD)
                tABCD = np.random.exponential(1.0)
                tCDpr = str(thetaABCD*tABCD)
                tABpr = str(thetaAB*(tau2-tAB)+thetaABC*tau3+thetaABCD*(tCD+tABCD))
                print '[1]((3:'+tApr+',4:'+tBpr+'):'+tABpr+',(2:'+tCpr+',1:'+tDpr+'):'+tCDpr+');'
                count[28] = count[28] + 1

    else: # A and B do not coalesce in their ancestral interval 
        time1 = np.random.exponential(1.0/3.0)  # next coalescent time
        if time1<tau3: # there is a coalescent event in population ABC
            ran1 = np.random.uniform()
            if ran1<1.0/3.0: # A and B coalesce first in population ABC
                tAB = time1
                tApr = str(thetaA*tau1+thetaAB*tau2+thetaABC*tAB)
                tBpr = str(thetaB*tau1+thetaAB*tau2+thetaABC*tAB)
                time2 = np.random.exponential(1.0)
                if time2<tau3-tAB: # AB and C coalesce in population ABC
                    tABC = time2
                    tABpr = str(thetaABC*tABC)  # tABpr = str(thetaABC*(tAB+tABC))
                    tCpr = str(thetaC*(tau1+tau2)+thetaABC*(tAB+tABC))
                    tABCD = np.random.exponential(1.0)
                    tABCpr = str(thetaABC*(tau3-(tAB+tABC))+thetaABCD*tABCD)
                    tDpr = str(thetaD*(tau1+tau2+tau3)+thetaABCD*tABCD)
                    print '[1](((3:'+tApr+',4:'+tBpr+'):'+tABpr+',2:'+tCpr+'):'+tABCpr+',1:'+tDpr+');'
                    count[2] = count[2]+1
                                
                else: #AB and C do not coalesce in population ABC
                    ran2 = np.random.uniform()
                    if ran2<1.0/3.0: #AB and C coalesce first, above the root
                        tABC = np.random.exponential(1.0/3.0)
                        tABpr = str(thetaABC*(tau3-tAB)+thetaABCD*tABC)
                        tCpr = str(thetaC*(tau1+tau2+tau3)+thetaABCD*tABC)
                        tABCD = np.random.exponential(1.0)
                        tABCpr = str(thetaABCD*tABCD)
                        tDpr = str(thetaD*(tau1+tau2+tau3)+thetaABCD*(tABC+tABCD))
                        print '[1](((3:'+tApr+',4:'+tBpr+'):'+tABpr+',2:'+tCpr+'):'+tABCpr+',1:'+tDpr+');'
                        count[3] = count[3]+1

                    elif ran2<2.0/3.0: # AB and D coalesce first, above the root
                        tABD = np.random.exponential(1.0/3.0)
                        tABpr = str(thetaABC*(tau3-tAB)+thetaABCD*tABD)
                        tDpr = str(thetaD*(tau1+tau2+tau3)+thetaABCD*tABD)
                        tABCD = np.random.exponential(1.0)
                        tABDpr = str(thetaABCD*tABCD)
                        tCpr = str(thetaC*(tau1+tau2)+thetaABC*tau3+thetaABCD*(tABD+tABCD))
                        print '[1](((3:'+tApr+',4:'+tBpr+'):'+tABpr+',1:'+tDpr+'):'+tABDpr+',2:'+tCpr+');'
                        count[27] = count[27] + 1

                    else: # C and D coalesce first, above the root
                        tCD = np.random.exponential(1.0/3.0)
                        tCpr = str(thetaC*(tau1+tau2)+thetaABC*tau3+thetaABCD*tCD)
                        tDpr = str(thetaD*(tau1+tau2+tau3)+thetaABCD*tCD)
                        tABCD = np.random.exponential(1.0)
                        tCDpr = str(thetaABCD*tABCD)
                        tABpr = str(thetaABC*(tau3-tAB)+thetaABCD*(tCD+tABCD))
                        print '[1]((3:'+tApr+',4:'+tBpr+'):'+tABpr+',(1:'+tDpr+',2:'+tCpr+'):'+tCDpr+');'
                        count[33] = count[33] + 1
                                   

            elif ran1<2.0/3.0: # A and C coalesce first in population ABC
                tAC = time1
                tApr = str(thetaA*tau1+thetaAB*tau2+thetaABC*tAC)
                tCpr = str(thetaC*(tau1+tau2)+thetaABC*tAC)
                time2 = np.random.exponential(1.0)
                if time2<tau3-tAC:  # AC and B coalesce in population ABC
                    tABC = time2
                    tACpr = str(thetaABC*tABC)
                    tBpr = str(thetaB*tau1+thetaAB*tau2+thetaABC*(tAC+tABC))
                    tABCD = np.random.exponential(1.0)
                    tABCpr = str(thetaABC*(tau3-(tAC+tABC))+thetaABCD*tABCD)
                    tDpr = str(thetaD*(tau1+tau2+tau3)+thetaABCD*tABCD)
                    print '[1](((3:'+tApr+',2:'+tCpr+'):'+tACpr+',4:'+tBpr+'):'+tABCpr+',1:'+tDpr+');'
                    count[5] = count[5]+1

                else: # AC and B do not coalesce in population ABC
                    ran1 = np.random.uniform()
                    if ran1<1.0/3.0: # AC and B coalesce first, above the root
                        tABC = np.random.exponential(1.0/3.0)
                        tACpr = str(thetaABC*(tau3-tAC)+ thetaABCD*tABC)
                        tBpr = str(thetaB*tau1+thetaAB*tau2+thetaABC*tau3+thetaABCD*tABC)
                        tABCD = np.random.exponential(1.0)
                        tABCpr = str(thetaABCD*tABCD)
                        tDpr = str(thetaD*(tau1+tau2+tau3)+thetaABCD*(tABC+tABCD))
                        print '[1](((3:'+tApr+',2:'+tCpr+'):'+tACpr+',4:'+tBpr+'):'+tABCpr+',1:'+tDpr+');'
                        count[6] = count[6]+1

                    elif ran1<2.0/3.0: # AC and D coalesce first, above the root
                        tACD = np.random.exponential(1.0/3.0)
                        tACpr = str(thetaABC*(tau3-tAC)+ thetaABCD*tACD)
                        tDpr = str(thetaD*(tau1+tau2+tau3)+thetaABCD*tACD)
                        tABCD = np.random.exponential(1.0)
                        tACDpr = str(thetaABCD*tABCD) # tACDpr = str(thetaABCD*tACD)
                        tBpr = str(thetaB*tau1+thetaAB*tau2+thetaABC*tau3+thetaABCD*(tACD+tABCD))
                        print '[1](((3:'+tApr+',2:'+tCpr+'):'+tACpr+',1:'+tDpr+'):'+tACDpr+',4:'+tBpr+');'
                        count[29] = count[29] + 1

                    else:  # B and D coalesce first, above the root
                        tBD = np.random.exponential(1.0/3.0)
                        tBpr = str(thetaB*tau1+thetaAB*tau2+thetaABC*tau3+thetaABCD*tBD)
                        tDpr = str(thetaD*(tau1+tau2+tau3)+thetaABCD*tBD)
                        tABCD = np.random.exponential(1.0)
                        tBDpr = str(thetaABCD*tABCD) # tBDpr = str(thetaABCD*tBD)
                        tACpr = str(thetaABC*(tau3-tAC)+thetaABCD*(tBD+tABCD))
                        print '[1]((3:'+tApr+',2:'+tCpr+'):'+tACpr+',(1:'+tDpr+',4:'+tBpr+'):'+tBDpr+');'
                        count[30] = count[30] + 1
                        
                    
            else: # B and C coalesce first in population ABC
                tBC = time1
                tBpr = str(thetaB*tau1+thetaAB*tau2+thetaABC*tBC)
                tCpr = str(thetaC*(tau1+tau2)+thetaABC*tBC)
                time2 = np.random.exponential(1.0)
                if time2<tau3-tBC: # BC and A coalesce in population ABC
                    tABC = time2
                    tBCpr = str(thetaABC*tABC)
                    tApr = str(thetaA*tau1+thetaAB*tau2+thetaABC*(tBC+tABC))
                    tABCD = np.random.exponential(1.0)
                    tABCpr = str(thetaABC*(tau3-(tBC+tABC))+thetaABCD*tABCD)
                    tDpr = str(thetaD*(tau1+tau2+tau3)+thetaABCD*tABCD)
                    print '[1](((4:'+tBpr+',2:'+tCpr+'):'+tBCpr+',3:'+tApr+'):'+tABCpr+',1:'+tDpr+');'
                    count[10] = count[10] + 1

                else: # BC and A do not coalesce in population ABC
                    ran1 = np.random.uniform()
                    if ran1<1.0/3.0: # BC and A coalesce first, above the root
                        tABC = np.random.exponential(1.0/3.0)
                        tBCpr  = str(thetaABC*(tau3-tBC)+thetaABCD*tABC)
                        tApr = str(thetaA*tau1+thetaAB*tau2+thetaABC*tau3+thetaABCD*tABC)
                        tABCD = np.random.exponential(1.0)
                        tABCpr = str(thetaABCD*tABCD)
                        tDpr = str(thetaD*(tau1+tau2+tau3)+thetaABCD*(tABC+tABCD))
                        print '[1](((4:'+tBpr+',2:'+tCpr+'):'+tBCpr+',3:'+tApr+'):'+tABCpr+',1:'+tDpr+');'
                        count[11] = count[11] + 1

                    elif ran1<2.0/3.0: # BC and D coalesce first, above the root
                        tBCD = np.random.exponential(1.0/3.0)
                        tBCpr = str(thetaABC*(tau3-tBC)+thetaABCD*tBCD)
                        tDpr = str(thetaD*(tau1+tau2+tau3)+thetaABCD*tBCD)
                        tABCD = np.random.exponential(1.0)
                        tBCDpr = str(thetaABCD*tABCD)
                        tApr = str(thetaA*tau1+thetaAB*tau2+thetaABC*tau3+thetaABCD*(tBCD+tABCD))
                        print '[1](((4:'+tBpr+',2:'+tCpr+'):'+tBCpr+',1:'+tDpr+'):'+tBCDpr+',3:'+tApr+');'
                        count[31] = count[31] + 1

                    else: # A and D coalesce first, above the root
                        tAD = np.random.exponential(1.0/3.0)
                        tApr = str(thetaA*tau1+thetaAB*tau2+thetaABC*tau3+thetaABCD*tAD)
                        tDpr = str(thetaD*(tau1+tau2+tau3)+thetaABCD*tAD)
                        tABCD = np.random.exponential(1.0)
                        tADpr = str(thetaABCD*tABCD)
                        tBCpr = str(thetaABC*(tau3-tBC)+thetaABCD*(tAD+tABCD))
                        print '[1]((4:'+tBpr+',2:'+tCpr+'):'+tBCpr+',(1:'+tDpr+',3:'+tApr+'):'+tADpr+');'
                        count[32] = count[32] + 1
                    

        else: # all coalescences above the root
            ran1 = np.random.uniform()
            if ran1<1.0/6.0: # A and B coalesce first
                tAB = np.random.exponential(1.0/6.0)
                tApr = str(thetaA*tau1+thetaAB*tau2+thetaABC*tau3+thetaABCD*tAB)
                tBpr = str(thetaB*tau1+thetaAB*tau2+thetaABC*tau3+thetaABCD*tAB)
                ran2 = np.random.uniform()
    
                if ran2<1.0/3.0:  # AB and C coalesce next
                    tABC = np.random.exponential(1.0/3.0)
                    tABpr = str(thetaABCD*tABC)
                    tCpr = str(thetaC*(tau1+tau2)+thetaABC*tau3+thetaABCD*(tAB+tABC))
                    tABCD = np.random.exponential(1.0)
                    tABCpr = str(thetaABCD*tABCD)
                    tDpr = str(thetaD*(tau1+tau2+tau3)+thetaABCD*(tAB+tABC+tABCD))
                    print '[1](((3:'+tApr+',4:'+tBpr+'):'+tABpr+',2:'+tCpr+'):'+tABCpr+',1:'+tDpr+');'
                    count[4] = count[4] + 1

                elif ran2<2.0/3.0: # AB and D coalesce next
                    tABD = np.random.exponential(1.0/3.0)
                    tABpr = str(thetaABCD*tABD)
                    tDpr = str(thetaD*(tau1+tau2+tau3)+thetaABCD*tABD)
                    tABCD = np.random.exponential(1.0)
                    tABDpr = str(thetaABCD*tABCD) # tABDpr = str(thetaABCD*tABD)
                    tCpr = str(thetaC*(tau1+tau2)+thetaABC*tau3+thetaABCD*(tAB+tABD+tABCD))
                    print '[1](((3:'+tApr+',4:'+tBpr+'):'+tABpr+',1:'+tDpr+'):'+tABDpr+',2:'+tCpr+');'
                    count[8] = count[8] + 1
                    
                else: # C and D coalesce next
                    tCD = np.random.exponential(1.0/3.0)
                    tCpr = str(thetaC*(tau1+tau2)+thetaABC*tau3+thetaABCD*tCD)
                    tDpr = str(thetaD*(tau1+tau2+tau3)+thetaABCD*tCD)
                    tABCD = np.random.exponential(1.0)
                    tABpr = str(thetaABCD*(tCD+tABCD))
                    tCDpr = str(thetaABCD*tABCD) # tCDpr = str(thetaABCD*tCD)
                    print '[1]((3:'+tApr+',4:'+tBpr+'):'+tABpr+',(2:'+tCpr+',1:'+tDpr+'):'+tCDpr+');'
                    count[9] = count[9] + 1

            elif ran1<2.0/6.0: # A and C coalesce first
                tAC = np.random.exponential(1.0/6.0)
                tApr = str(thetaA*tau1+thetaAB*tau2+thetaABC*tau3+thetaABCD*tAC)
                tCpr = str(thetaC*(tau1+tau2)+thetaABC*tau3+thetaABCD*tAC)
                ran2 = np.random.uniform()
                if ran2<1.0/3.0:  # AC and B coalesce next
                    tABC = np.random.exponential(1.0/3.0)
                    tACpr = str(thetaABCD*tABC)
                    tBpr = str(thetaB*tau1+thetaAB*tau2+thetaABC*tau3+thetaABCD*(tAB+tABC))
                    tABCD = np.random.exponential(1.0)
                    tABCpr = str(thetaABCD*tABCD)
                    tDpr = str(thetaD*(tau1+tau2+tau3)+thetaABCD*(tAB+tABC+tABCD))
                    print '[1](((3:'+tApr+',2:'+tCpr+'):'+tACpr+',4:'+tBpr+'):'+tABCpr+',1:'+tDpr+');'
                    count[7] = count[7] + 1

                elif ran2<2.0/3.0: # AC and D coalesce next
                    tACD = np.random.exponential(1.0/3.0)
                    tACpr = str(thetaABCD*tACD)
                    tDpr = str(thetaD*(tau1+tau2+tau3)+thetaABCD*(tAC+tACD))
                    tABCD = np.random.exponential(1.0)
                    tACDpr = str(thetaABCD*tABCD)
                    tBpr = str(thetaB*tau1+thetaAB*tau2+thetaABC*tau3+thetaABCD*(tAC+tACD+tABCD))
                    print '[1](((3:'+tApr+',2:'+tCpr+'):'+tACpr+',1:'+tDpr+'):'+tACDpr+',4:'+tBpr+');'
                    count[15] = count[15] + 1

                else: # B and D coalesce next
                    tBD = np.random.exponential(1.0/3.0)
                    tBpr = str(thetaB*tau1+thetaAB*tau2+thetaABC*tau3+thetaABCD*(tAC+tBD))
                    tDpr = str(thetaD*(tau1+tau2+tau3)+thetaABCD*(tAC+tBD))
                    tABCD = np.random.exponential(1.0)
                    tACpr = str(thetaABCD*(tBD+tABCD))
                    tBDpr = str(thetaABCD*tABCD)
                    print '[1]((3:'+tApr+',2:'+tCpr+'):'+tACpr+',(4:'+tBpr+',1:'+tDpr+'):'+tBDpr+');'
                    count[16] = count[16] + 1

            elif ran1<3.0/6.0: # B and C coalesce first
                tBC = np.random.exponential(1.0/6.0)
                tBpr = str(thetaB*tau1+thetaAB*tau2+thetaABC*tau3+thetaABCD*tBC)
                tCpr = str(thetaC*(tau1+tau2)+thetaABC*tau3+thetaABCD*tBC)
                ran2 = np.random.uniform()
                if ran2<1.0/3.0: # BC and A coalesce next
                    tABC = np.random.exponential(1.0/3.0)
                    tBCpr = str(thetaABCD*tABC)
                    tApr = str(thetaA*tau1+thetaAB*tau2+thetaABC*tau3+thetaABCD*(tAB+tABC))
                    tABCD = np.random.exponential(1.0)
                    tABCpr = str(thetaABCD*tABCD)
                    tDpr = str(thetaD*(tau1+tau2+tau3)+thetaABCD*(tAB+tABC+tABCD))
                    print '[1](((4:'+tBpr+',2:'+tCpr+'):'+tBCpr+',3:'+tApr+'):'+tABCpr+',1:'+tDpr+');'
                    count[12] = count[12] + 1

                elif ran2<2.0/3.0: # BC and D coalesce next
                    tBCD = np.random.exponential(1.0/3.0)
                    tBCpr = str(thetaABCD*tBCD)
                    tDpr = str(thetaD*(tau1+tau2+tau3)+thetaABCD*(tBC+tBCD))
                    tABCD = np.random.exponential(1.0)
                    tBCDpr = str(thetaABCD*tABCD)
                    tApr = str(thetaA*tau1+thetaAB*tau2+thetaABC*tau3+thetaABCD*(tBC+tBCD+tABCD))
                    print '[1](((4:'+tBpr+',2:'+tCpr+'):'+tBCpr+',1:'+tDpr+'):'+tBCDpr+',3:'+tApr+');'
                    count[13] = count[13] + 1

                else: # A and D coalesce next
                    tAD = np.random.exponential(1.0/3.0)
                    tApr = str(thetaA*tau1+thetaAB*tau2+thetaABC*tau3+thetaABCD*(tBC+tAD))
                    tDpr = str(thetaD*(tau1+tau2+tau3)+thetaABCD*(tBC+tAD))
                    tABCD = np.random.exponential(1.0)
                    tADpr = str(thetaABCD*tABCD)
                    tBCpr = str(thetaABCD*(tAD+tABCD))
                    print '[1]((4:'+tBpr+',2:'+tCpr+'):'+tBCpr+',(3:'+tApr+',1:'+tDpr+'):'+tADpr+');'
                    count[14] = count[14] + 1

            elif ran1<4.0/6.0: # A and D coalesce first
                tAD = np.random.exponential(1.0/6.0)
                tApr = str(thetaA*tau1+thetaAB*tau2+thetaABC*tau3+thetaABCD*tAD)
                tDpr = str(thetaD*(tau1+tau2+tau3)+thetaABCD*tAD)
                ran2 = np.random.uniform()
                if ran2<1.0/3.0:  # AD and B coalesce next
                    tABD = np.random.exponential(1.0/3.0)
                    tADpr = str(thetaABCD*tABD)
                    tBpr = str(thetaB*tau1+thetaAB*tau2+thetaABC*tau3+thetaABCD*(tAD+tABD))
                    tABCD = np.random.exponential(1.0)
                    tABDpr = str(thetaABCD*tABCD)
                    tCpr = str(thetaC*(tau1+tau2)+thetaABC*tau3+thetaABCD*(tAD+tABD+tABCD))
                    print '[1](((3:'+tApr+',1:'+tDpr+'):'+tADpr+',4:'+tBpr+'):'+tABDpr+',2:'+tCpr+');'
                    count[17] = count[17] + 1

                elif ran2<2.0/3.0: # AD and C coalesce next
                    tACD = np.random.exponential(1.0/3.0)
                    tADpr = str(thetaABCD*tACD)
                    tCpr = str(thetaC*(tau1+tau2)+thetaABC*tau3+thetaABCD*(tAD+tACD))
                    tABCD = np.random.exponential(1.0)
                    tACDpr = str(thetaABCD*tABCD) # tACDpr = str(thetaABCD*tACD)
                    tBpr = str(thetaB*tau1+thetaAB*tau2+thetaABC*tau3+thetaABCD*(tAD+tACD+tABCD))
                    print '[1](((3:'+tApr+',1:'+tDpr+'):'+tADpr+',2:'+tCpr+'):'+tACDpr+',4:'+tBpr+');'
                    count[18] = count[18] + 1

                else: # B and C coalesce next
                    tBC = np.random.exponential(1.0/3.0)
                    tBpr = str(thetaB*tau1+thetaAB*tau2+thetaABC*tau3+thetaABCD*tBC)
                    tCpr = str(thetaC*(tau1+tau2)+thetaABC*tau3+thetaABCD*tBC)
                    tABCD = np.random.exponential(1.0)
                    tBCpr = str(thetaABCD*tABCD)
                    tADpr = str(thetaABCD*(tBC+tABCD))
                    print '[1]((3:'+tApr+',1:'+tDpr+'):'+tADpr+',(2:'+tCpr+',4:'+tBpr+'):'+tBCpr+');'
                    count[19] = count[19] + 1

            elif ran1<5.0/6.0: # B and D coalesce first
                tBD = np.random.exponential(1.0/6.0)
                tBpr = str(thetaB*tau1+thetaAB*tau2+thetaABC*tau3+thetaABCD*tBD)
                tDpr = str(thetaD*(tau1+tau2+tau3)+thetaABCD*tBD)
                ran2 = np.random.uniform()
                if ran2<1.0/3.0: # BD and A coalesce next
                    tABD = np.random.exponential(1.0/6.0)
                    tBDpr = str(thetaABCD*tABD)
                    tApr = str(thetaA*tau1+thetaAB*tau2+thetaABC*tau3+thetaABCD*(tBD+tABD))
                    tABCD = np.random.exponential(1.0)
                    tABDpr = str(thetaABCD*tABCD)
                    tCpr = str(thetaC*(tau1+tau2)+thetaABC*tau3+thetaABCD*(tBD+tABD+tABCD))
                    print '[1](((4:'+tBpr+',1:'+tDpr+'):'+tBDpr+',3:'+tApr+'):'+tABDpr+',2:'+tCpr+');'
                    count[20] = count[20] + 1

                elif ran2<2.0/3.0: # BD and C coalesce next
                    tBCD = np.random.exponential(1.0/3.0)
                    tBDpr = str(thetaABCD*tBCD)
                    tCpr = str(thetaC*(tau1+tau2)+thetaABC*tau3+thetaABCD*(tBD+tBCD))
                    tABCD = np.random.exponential(1.0)
                    tBCDpr = str(thetaABCD*tABCD)
                    tApr = str(thetaA*tau1+thetaAB*tau2+thetaABC*tau3+thetaABCD*(tBD+tBCD+tABCD))
                    print '[1](((4:'+tBpr+',1:'+tDpr+'):'+tBDpr+',2:'+tCpr+'):'+tBCDpr+',3:'+tApr+');'
                    count[21] = count[21] + 1

                else: # A and C coalesce next
                    tAC = np.random.exponential(1.0/3.0)
                    tApr = str(thetaA*tau1+thetaAB*tau2+thetaABC*tau3+thetaABCD*(tBD+tAC))
                    tCpr = str(thetaC*(tau1+tau2)+thetaABC*tau3+thetaABCD*(tBD+tAC))
                    tABCD = np.random.exponential(1.0)
                    tACpr = str(thetaABCD*tABCD)
                    tBDpr = str(thetaABCD*(tAC+tABCD))
                    print '[1]((4:'+tBpr+',1:'+tDpr+'):'+tBDpr+',(2:'+tCpr+',3:'+tApr+'):'+tACpr+');'
                    count[22] = count[22] + 1

            else: # C and D coalesce first
                tCD = np.random.exponential(1.0/6.0)
                tCpr = str(thetaC*(tau1+tau2)+thetaABC*tau3+thetaABCD*tCD)
                tDpr = str(thetaD*(tau1+tau2+tau3)+thetaABCD*tCD)
                ran2 = np.random.uniform()
                if ran2<1.0/3.0:  # CD and A coalesce next
                    tACD = np.random.exponential(1.0/3.0)
                    tCDpr = str(thetaABCD*tACD)
                    tApr = str(thetaA*tau1+thetaAB*tau2+thetaABC*tau3+thetaABCD*(tCD+tACD))
                    tABCD = np.random.exponential(1.0)
                    tACDpr = str(thetaABCD*tABCD)
                    tBpr = str(thetaB*tau1+thetaAB*tau2+thetaABC*tau3+thetaABCD*(tCD+tACD+tABCD))
                    print '[1](((2:'+tCpr+',1:'+tDpr+'):'+tCDpr+',3:'+tApr+'):'+tACDpr+',4:'+tBpr+');'
                    count[23] = count[23] + 1

                elif ran2<2.0/3.0: # CD and B coalesce next
                    tBCD = np.random.exponential(1.0/3.0)
                    tCDpr = str(thetaABCD*tBCD)
                    tBpr = str(thetaB*tau1+thetaAB*tau2+thetaABC*tau3+thetaABCD*tBCD)
                    tABCD = np.random.exponential(1.0)
                    tBCDpr = str(thetaABCD*tABCD)
                    tApr = str(thetaA*tau1+thetaAB*tau2+thetaABC*tau3+thetaABCD*(tCD+tBCD+tABCD))
                    print '[1](((2:'+tCpr+',1:'+tDpr+'):'+tCDpr+',4:'+tBpr+'):'+tBCDpr+',3:'+tApr+');'
                    count[24] = count[24] + 1

                else: # A and B coalesce next
                    tAB =np.random.exponential(1.0/3.0)
                    tApr = str(thetaA*tau1+thetaAB*tau2+thetaABC*tau3+thetaABCD*tAB)
                    tBpr = str(thetaB*tau1+thetaAB*tau2+thetaABC*tau3+thetaABCD*tAB)
                    tABCD = np.random.exponential(1.0)
                    tABpr = str(thetaABCD*tABCD)
                    tCDpr = str(thetaABCD*(tAB+tABCD))
                    print '[1]((2:'+tCpr+',1:'+tDpr+'):'+tCDpr+',(4:'+tBpr+',3:'+tApr+'):'+tABpr+');'
                    count[25] = count[25] + 1
                    
                    

                    

                    
                    
            


#print'\n\n'
#print 'History probabilities'
#print count/ntrees
#print '\n'
#print sum(count)/ntrees

#print '(((A,B),C),D)',(count[0]+count[1]+count[2]+count[3]+count[4])/ntrees
#print '(((A,C),B),D)',(count[5]+count[6]+count[7])/ntrees
#print '(((A,B),D),C)',(count[8]+count[26]+count[27])/ntrees
#print '((A,B),(D,C))',(count[9]+count[25]+count[28]+count[33])/ntrees
#print '(((B,C),A),D)',(count[10]+count[11]+count[12])/ntrees
#print '(((B,C),D),A)',(count[13]+count[31])/ntrees
#print '(((A,D),B),C)',(count[17])/ntrees
#print '((B,C),(A,D))',(count[14]+count[19]+count[32])/ntrees
#print '(((A,C),D),B)',(count[15]+count[29])/ntrees
#print '((B,D),(A,C))',(count[16]+count[22]+count[30])/ntrees
#print '(((A,D),C),B)',(count[18])/ntrees
#print '(((B,D),A),C)',(count[20])/ntrees
#print '(((B,D),C),A)',(count[21])/ntrees
#print '(((C,D),A),B)',(count[23])/ntrees
#print '(((C,D),B),A)',(count[24])/ntrees
