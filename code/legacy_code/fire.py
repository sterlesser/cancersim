import scipy as sp

def FIRE(x0,fprime,fmax=0.005,
            Nmin=5.,finc=1.1,fdec=0.5,alphastart=0.1,fa=0.99,deltatmax=10.,
            maxsteps = 10**5):

    Nmin,finc,fdec,alphastart,fa,deltatmax=(5.,1.1,0.5,0.1,0.99,10.)

    alpha = alphastart
    deltat = 0.1

    pos = x0.copy()
    v = sp.zeros_like(pos)

    steps_since_negative = 0

    def norm(vec):
        return sp.sqrt(sp.sum(vec**2,1))
    def unitize(vec):
        return ((vec.T)/norm(vec)).T

    forces = fprime(pos)

    step_num = 0

    while max(norm(forces)) > fmax and step_num < maxsteps:
        forces = fprime(pos)
        power = sp.vdot(forces,v)

        print "Step: {}, max_force: {}, power: {}".format(step_num,max(norm(forces)), power)

        v = (1.0 - alpha)*v + alpha*(norm(v)*unitize(forces).T).T

        if power>0.:
            if steps_since_negative > Nmin:
                deltat = min(deltat * finc, deltatmax)
                alpha = alpha*fa
            steps_since_negative += 1

        else:
            steps_since_negative = 0

            deltat = deltat * fdec
            v *= 0.
            alpha = alphastart

        v += forces*deltat
        pos += v*deltat
        step_num += 1
    
    return pos

    