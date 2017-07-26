#!/usr/bin/env python
# ! ## File: dynSIS.py
# ! Module: use networkx graphs!
# ! ## See README.md for more information and use
# !-----------------------------------------------------------------------------
# ! SIS epidemic model algorithm based on the article
# !           Computer Physics Communications 219C (2017) pp. 303-312
# !           "Optimized Gillespie algorithms for the simulation of 
# !            Markovian epidemic processes on large and heterogeneous networks"
# ! Copyright (C) 2017 Wesley Cota, Silvio C. Ferreira
# ! 
# ! Please cite the above cited paper (available at <http://dx.doi.org/10.1016/j.cpc.2017.06.007> ) 
# ! as reference to our code.
# ! 
# !    This program is free software: you can redistribute it and/or modify
# !    it under the terms of the GNU General Public License as published by
# !    the Free Software Foundation, either version 3 of the License, or
# !    (at your option) any later version.
# !
# !    This program is distributed in the hope that it will be useful,
# !    but WITHOUT ANY WARRANTY; without even the implied warranty of
# !    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# !    GNU General Public License for more details.
# !
# !    You should have received a copy of the GNU General Public License
# !    along with this program.  If not, see <http://www.gnu.org/licenses/>.
# !-----------------------------------------------------------------------------
# ! Author    : Wesley Cota
# ! Email     : wesley.cota@ufv.br
# ! Date      : 27 Mar 2017
# ! Version   : 1.0
# !-----------------------------------------------------------------------------
# ! See README.md for more details
# ! This code is available at <https://github.com/wcota/dynSIS-networkx>
# ! For performance, see <https://github.com/wcota/dynSIS> (Fortran implementation)
# ! For pure Python, see <https://github.com/wcota/dynSIS-py>

import numpy as np
from math import log

print(  '################################################################################',
        '######### Optimized Gillespie algorithms for the simulation of Markovian  ######',
        '####### epidemic processes on large and heterogeneous networks: SIS-OGA. #######',
        '##============ Copyright (C) 2017 Wesley Cota, Silvio C. Ferreira ============##',
        '##===== Paper available at <http://dx.doi.org/10.1016/j.cpc.2017.06.007> =====##',
        '##======= The codes are available at <https://github.com/wcota/dynSIS> =======##',
        '##======== Please cite the above cited paper as reference to our code ========##',
        '##=== This code is under GNU General Public License. Please see README.md. ===##',
        '################################################################################',
        '',
        sep='\n')

def dyn_run(nw, fnOutput, dynp_sam, dynp_lb, dynp_tmax, dynp_pINI):
    net_N = nw.number_of_nodes()
    net_kmax = max(nw.degree().values())        # Used in the rejection probability
    avg_rho = np.zeros(dynp_tmax, np.float64)   # Average for rho at times t, averaged
    avg_t = np.zeros(dynp_tmax, np.float64)
    avg_sam = np.zeros(dynp_tmax, np.int)       # number of samples for each time t
    avg_samSurv = np.zeros(dynp_tmax, np.int)   # and of survivng ones

    dyn_VI = [None]*net_N                       # list V^I. Any node type is allowed
    dyn_sig = { i : 0 for i in nw.nodes()}      # sigma
    
    print('\nOk! Doing the samples...')
    
    dyn_dt_pos_max = 0
    for sam in range(1,dynp_sam+1):
        print('\nSample #', sam)
        
        print('|| (random) Initial condition...')
        
        dyn_sig = dict.fromkeys(dyn_sig, 0)
        dyn_VI = [None]*net_N
        dyn_NI = 0
        dyn_Nk = 0
        
        # Sort vertices and apply the initial condition
        for ver in np.random.permutation(nw.nodes()):
            dyn_VI[dyn_NI] = ver
            dyn_NI += 1
            dyn_sig[ver]= 1
            dyn_Nk += nw.degree(ver)
            if dyn_NI == int(net_N*dynp_pINI):
                break
    
        # Run dynamics
        dyn_t = 0
        dyn_dt = 0.0
        dyn_dt_pos = 1
    
        print('|| Running dynamics...')
        while dyn_t <= dynp_tmax and dyn_NI > 0:
            # SIS-OGA ALGORITHM
        
            # Calculate the total rate
            dyn_R = (dyn_NI + 1.0*dynp_lb * dyn_Nk)
        
            # Select the time step
            rnd = max(np.random.uniform(),1e-12) # Avoid u = 0
            dyn_dt = -log(rnd) / dyn_R
        
            # Update the time
            dyn_t += dyn_dt
        
            # Probability m to heal
            dyn_m = 1.0*dyn_NI / dyn_R
        
            if np.random.uniform() < dyn_m: # Select a random occupied vertex and heal.
                pos_inf = np.random.randint(0,dyn_NI)
                ver = dyn_VI[pos_inf]
                
                # Then, heal it
                dyn_sig[ver] = 0
                dyn_Nk -= nw.degree(ver)
                dyn_NI -= 1
                dyn_VI[pos_inf] = dyn_VI[dyn_NI]
            else: # If not, try to infect: w = 1 - m
                # Select the infected vertex i with prob. proportional to k_i
                while True:
                    pos_inf = np.random.randint(0,dyn_NI)
                    ver = dyn_VI[pos_inf]
                    if np.random.uniform() < 1.0*nw.degree(ver) / (1.0*net_kmax):
                        break
            
                # Select one of its neighbors
                ver = np.random.choice(nw.neighbors(ver))
            
                if dyn_sig[ver] == 0: # if not a phantom process, infect
                    dyn_sig[ver] = 1
                    dyn_Nk += nw.degree(ver)
                    dyn_VI[dyn_NI] = ver    # Add one element to list
                    dyn_NI += 1             # Increase by 1 the list
    
                # Try to save the dynamics by time unit
                while (dyn_t >= dyn_dt_pos): # Save data
                    avg_rho[dyn_dt_pos - 1] += 1.0*dyn_NI/net_N
                    avg_t[dyn_dt_pos - 1] += dyn_t
                    avg_sam[dyn_dt_pos - 1] += 1
                    if dyn_NI != 0: 
                        avg_samSurv[dyn_dt_pos - 1] += 1
                        dyn_dt_pos_max = max(dyn_dt_pos,dyn_dt_pos_max) # The maximum t with non-null rho
                    dyn_dt_pos += 1
                
                # if a absorbing state is reached, exit
    
        # Write output file
        flOutput = open(fnOutput, 'wt')
        print(  '## ***** Algorithm used: Optimized Gillespie Algorithm for SIS (SIS-OGA, NetworkX) *****',
                '#@ Number of nodes: '+str(net_N),
                '#@ Number of edges: '+str(2*nw.number_of_edges()),
                '#@ Samples: '+str(dynp_sam),
                '#! Infection rate (lambda): '+str(dynp_lb),
                '#! Maximum time steps: '+str(dynp_tmax),
                '#! Fraction of infected vertices (initial condition): '+str(dynp_pINI),
                sep='\n',
                file=flOutput)
            
        for dt_pos in range(0,dyn_dt_pos_max):
            print(1.0*avg_t[dt_pos]/avg_sam[dt_pos], 1.0*avg_rho[dt_pos]/(1.0*sam),
                    file=flOutput)
        # If you use /avg_samSurv[dt_pos] instead of /(1.0*sam) to write avg_rho (2nd column), you have 
        # QS analysis :)
                
        flOutput.close()
