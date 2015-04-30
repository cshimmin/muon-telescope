#!/usr/bin/env python

import numpy as np
from scipy import stats

def zenith_distribution_discrete(n=2):
    ''' Return a discrete distribution representing a
    zenith distribution of cos(x)^n
    '''
    # more points = finer angular scales resolved
    N_points = 1000 
    # generate the bins over the range [0,pi/2]
    X = np.linspace(0,np.pi/2.,N_points)
    A = np.arange(N_points)
    # sample the cos^n distribution and normalize the bins
    pk = np.cos(X)**n
    pk /= np.sum(pk)
    class rv(stats.rv_discrete):
        def rvs(self, *args, **kwargs):
            # for some reason scipy's rv_discrete can only handle
            # integers, so we have to translate the bins to a
            # value between [0, pi/2] manually.
            result = stats.rv_discrete.rvs(self, *args, **kwargs)
            return X[result]
    return rv(values=(A, pk))

def get_track_angle(zen_dist):
    ''' get a random (theta, phi) angle air defining a track direction,
    drawn from the given zenith distribution. '''
    theta = zen_dist.rvs()
    phi = np.random.uniform(0,2.*np.pi)
    return theta, phi

def get_unit_vector(theta, phi):
    ''' get a unit 3-vector corresponding to the given (phi, theta) angles '''
    sth = np.sin(theta)
    cth = np.cos(theta)
    return np.array([np.cos(phi)*sth, np.sin(phi)*sth, cth])

if __name__ == "__main__":
    scint_size = 3.16 # cm
    scint_spacing = 5 # cm
    n_paddles = 3

    n_trials = 100000

    zen_dist = zenith_distribution_discrete(n=2)

    hit_record = [] # keep a record of which layer(s) were hit in each event
    for itrial in xrange(n_trials):
        if itrial%(n_trials/10)==0:
            print "sampling %d/%d (%.1f%%)" % (itrial, n_trials, 100.*itrial/n_trials)
        # randomly pick a paddle layer that is hit by a cosmic ray
        ipad0 = np.random.randint(n_paddles)

        # pick a random location on that paddle
        x0, y0 = np.random.uniform(-scint_size/2., scint_size/2., 2)

        # next generate a track passing through this location
        theta0, phi0 = get_track_angle(zen_dist)
        v0 = get_unit_vector(theta0, phi0)

        # now check which of the paddles were hit by this track

        hit_code = 0 # bits representing which layer(s) were hit
        for ipad in xrange(n_paddles):
            if ipad == ipad0:
                # by definition this paddle was hit!
                hit_code |= 2**ipad
                continue
            
            # project the track to this layer
            vertical_distance = (ipad - ipad0) * scint_spacing
            v = v0 * vertical_distance / v0[2]

            # shift by the original x,y offset
            v[:2] += [x0, y0]

            # and check whether it is in the bounds of this layer
            layer_hit = np.all(np.abs(v[:2]) < scint_size/2.)
            if layer_hit:
                hit_code |= 2**ipad

        hit_record.append(hit_code)

    # print stats about the acceptance at each layer (assuming 3 layer problem)
    hit_record = np.array(hit_record)
    abc = np.sum( hit_record == 7 )
    ab = np.sum( (hit_record & 3) == 3 )
    bc = np.sum( (hit_record & 6) == 6 )
    ac = np.sum( (hit_record & 5) == 5 )

    print "N_abc =", abc
    print "N_abc/Nab =", (1.*abc/ab)
    print "N_abc/Nbc =", (1.*abc/bc)
    print "N_abc/Nac =", (1.*abc/ac)
