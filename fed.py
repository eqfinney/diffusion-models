import numpy as np

from kappa_lookup import kappa_lookup

def is_prime(n):
    s = np.sqrt(n+1.0)
    for i in xrange(2, int(s+1)):
        if n % i == 0:
            return False
    return True


def _fed_tau_internal(num_time_steps, scale, stability, reorder_flag):
    tau = np.zeros(num_time_steps)
    
    # saving future computations
    denom = 1./(4.*num_time_steps + 2.)
    num = scale * stability / 2.
    
    # setting up tau vector
    for k in xrange(num_time_steps):
        h = np.cos(np.pi * ( (2.*k + 1.) * denom ) )
        tau[k] = num / (h*h)
    
    if reorder_flag == True:
        permuted_indices = get_permutation(num_time_steps)
        tau = tau[permuted_indices]
            
    return num_time_steps, tau


def get_permutation(num_time_steps):
    if num_time_steps < len(kappa_lookup): # FED_MAXKAPPA in the C code
        kappa = kappa_lookup[num_time_steps]
    else: 
        kappa = int(num_time_steps/4)

    # find the prime on which to permute
    prime = num_time_steps + 1
    while not is_prime(prime):
        prime += 1

    # perform the permutation!
    indices = []
    k = 0
    for l in xrange(num_time_steps):
        index = (k+1)*kappa % prime - 1.
        while index >= num_time_steps:
            k += 1
            index = (k+1)*kappa % prime - 1.
        indices.append(int(index))
        k+=1
        
    return indices


def fed_tau_by_process_time(stop_time, num_cycles, stability, reorder_flag):
    cycle_time = stop_time/num_cycles
    N, tau = fed_tau_by_cycle_time(cycle_time, stability, reorder_flag)
    return N, tau


def fed_tau_by_cycle_time(stop_time, stability, reorder_flag):
    num_time_steps = int(np.ceil(np.sqrt(3.*stop_time / stability + 0.25)
                                 -0.5 - 1.e-8) + 0.5)
    scale = 3.*stop_time / (stability*float(num_time_steps*(num_time_steps
                                                            +1)))
    N, tau = _fed_tau_internal(num_time_steps, scale, stability, 
                               reorder_flag)
    return N, tau


