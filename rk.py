import math

def rk4(function, *args, **kwargs):

    t = float(args[0])
    neq = len(args[1]) # y - state vector
    h = float(args[2])
    k1 = [None] * neq
    t1 = [None] * neq
    k2 = [None] * neq
    t2 = [None] * neq
    k3 = [None] * neq
    t3 = [None] * neq
    k4 = [None] * neq

    for i in range(0, neq):
        k1[i] = h * function(t, args[1], i)
    for i in range(0, neq):
        t1[i] = args[1][i] + 0.5 * k1[i]
    for i in range(0, neq):
        k2[i] = h * function(t + 0.5 * h, t1, i)
    for i in range(0, neq):
        t2[i] = args[1][i] + 0.5 * k2[i];
    for i in range(0, neq):
        k3[i] = h * function(t + 0.5 * h, t2, i)
    for i in range(0, neq):
        t3[i] = args[1][i] + k3[i]
    t += h
    for i in range(0, neq):
        k4[i] = h * function(t, t3, i);
    for i in range(0, neq):
        args[1][i] += 1.0/6.0 * ( k1[i] + 2.0 * k2[i]+ 2.0 * k3[i] + k4[i] );


def rk45(function, *args, **kwargs):

    t = float(args[0])
    neq = len(args[1]) # y - state vector
    h = float(args[2])
    tol = float(args[3])

    k1 = [None] * neq
    k2 = [None] * neq
    t2 = [None] * neq
    k3 = [None] * neq
    t3 = [None] * neq
    k4 = [None] * neq
    t4 = [None] * neq
    k5 = [None] * neq
    t5 = [None] * neq
    k6 = [None] * neq
    t6 = [None] * neq
    nv = [None] * neq
    z  = [None] * neq

    # Coefficients used to compute the independent variable argument of f

    c20  =   2.500000000000000e-01  #  1/4
    c30  =   3.750000000000000e-01  #  3/8
    c40  =   9.230769230769231e-01  #  12/13
    c50  =   1.000000000000000e+00  #  1
    c60  =   5.000000000000000e-01  #  1/2

    # Coefficients used to compute the dependent variable argument of f

    c21 =   2.500000000000000e-01  #  1/4
    c31 =   9.375000000000000e-02  #  3/32
    c32 =   2.812500000000000e-01  #  9/32
    c41 =   8.793809740555303e-01  #  1932/2197
    c42 =  -3.277196176604461e+00  # -7200/2197
    c43 =   3.320892125625853e+00  #  7296/2197
    c51 =   2.032407407407407e+00  #  439/216
    c52 =  -8.000000000000000e+00  # -8
    c53 =   7.173489278752436e+00  #  3680/513
    c54 =  -2.058966861598441e-01  # -845/4104
    c61 =  -2.962962962962963e-01  # -8/27
    c62 =   2.000000000000000e+00  #  2
    c63 =  -1.381676413255361e+00  # -3544/2565
    c64 =   4.529727095516569e-01  #  1859/4104
    c65 =  -2.750000000000000e-01  # -11/40

    # Coefficients used to compute 4th order RK estimate

    a1  =   1.157407407407407e-01  #  25/216
    a2  =   0.000000000000000e-00  #  0
    a3  =   5.489278752436647e-01  #  1408/2565
    a4  =   5.353313840155945e-01  #  2197/4104
    a5  =  -2.000000000000000e-01  # -1/5

    b1  =   1.185185185185185e-01  #  16.0/135.0
    b2  =   0.000000000000000e-00  #  0
    b3  =   5.189863547758284e-01  #  6656.0/12825.0
    b4  =   5.061314903420167e-01  #  28561.0/56430.0
    b5  =  -1.800000000000000e-01  # -9.0/50.0
    b6  =   3.636363636363636e-02  #  2.0/55.0

    for i in range(0, neq):
        k1[i] = h * function(t, args[1], i)

    for i in range(0, neq):
        t2[i] = args[1][i] + c20*k1[i]
    for i in range(0, neq):
        k2[i] = h * function(t + c20*h, t2, i)

    for i in range(0, neq):
        t3[i] = args[1][i] + c31*k1[i] + c32*k2[i]
    for i in range(0, neq):
        k3[i] = h * function(t + c30*h, t3, i)

    for i in range(0, neq):
        t4[i] = args[1][i] + c41*k1[i] + c42*k2[i] + c43*k3[i]
    for i in range(0, neq):
        k4[i] = h * function(t + c40*h , t4, i)

    for i in range(0, neq):
        t5[i] = args[1][i] + c51*k1[i] + c52*k2[i] + c53*k3[i] + c54*k4[i]
    for i in range(0, neq):
        k5[i] = h * function(t + c50*h, t5, i)

    for i in range(0, neq):
        t6[i] = args[1][i] + c61*k1[i] + c62*k2[i] + c63*k3[i] + c64*k4[i] + c65*k5[i]
    for i in range(0, neq):
        k6[i] = h * function(t + c60*h, t6, i)

    for i in range(0, neq):
        nv[i] = args[1][i] + a1*k1[i] + a2*k2[i] + a3*k3[i] + a4*k4[i] + a5*k5[i]

    for i in range(0, neq):
        z[i] = args[1][i] + b1*k1[i] + b2*k2[i] + b3*k3[i] + b4*k4[i] + b5*k5[i] + b6*k6[i]

    diff = 0
    for i in range(0, neq):
        diff +=  (nv[i] - z[i])*(nv[i] - z[i])
    diff = math.sqrt(diff)
    if diff==0:
        diff=1e-16

#    newstep = h * 0.84*math.pow (tol/diff, 0.25)
    newstep = h * min( max( 0.84 * ( tol * h / diff )**0.25, 0.1 ), 4.0 )
#    print (diff, tol, math.sqrt(nv[0]*nv[0] + nv[1]*nv[1] + nv[2]*nv[2] + nv[3]*nv[3]) )

    for i in range(0, neq):
        args[1][i] = nv[i]

    return newstep