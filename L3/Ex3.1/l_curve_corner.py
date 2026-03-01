def l_curve_corner(rho,eta,reg_param):
    # Parameter Estimation and Inverse Problems, 2nd edition, 2011
    # by R. Aster, B. Borchers, C. Thurber
    #
    # reg_corner, ireg_corner, kappa = l_curve_corner(rho,eta,reg_param)
    # returns l curve corner estimated using a maximum curvature (kappa) estimation 
    # in log-log space
    # rho is the misfit and eta is the model norm or seminorm
    #
    # INPUT
    #   rho       - misfit
    #   eta       - model norm or seminorm
    #   reg_param - the regularization parameter
    #
    # OUTPUT
    #   reg_corner  - the value of reg_param with maximum curvature
    #   ireg_corner - the index of the value in reg_param with maximum curvature
    #   kappa       - the curvature for each reg_param
    #
    #transform rho and eta into log-log space
    rho = np.array(rho, dtype = np.float64)
    eta = np.array(eta, dtype = np.float64)
    reg_param = np.array(reg_param, dtype = np.float64)
    
    x=np.log10(rho)
    y=np.log10(eta)

    # Triangular/circumscribed circle simple approximation to curvature 
    # (after Roger Stafford)

    # the series of points used for the triangle/circle
    x1 = x[0:-3]
    x2 = x[1:-2]
    x3 = x[2:-1]
    y1 = y[0:-3]
    y2 = y[1:-2]
    y3 = y[2:-1]

    # the side lengths for each triangle
    a = np.sqrt((x3-x2)**2+(y3-y2)**2)
    b = np.sqrt((x1-x3)**2+(y1-y3)**2)
    c = np.sqrt((x2-x1)**2+(y2-y1)**2)

    s=(a+b+c)/2 # semi-perimeter

    # the radius of each circle
    sa = s-a
    sb = s-b
    sc = s-c

    C_div = 4*np.sqrt((s*(sa)*(sb)*(sc)))
    R=(a*b*c)/C_div
    
    
    
    # The curvature for each estimate for each value which is
    # the reciprocal of its circumscribed radius. Since there aren't circles for 
    # the end points they have no curvature
    kappa = np.hstack((0, 1/R, 0)).T
    ireg_corner=np.nanargmax(abs(kappa))
    reg_corner=reg_param[ireg_corner]
    
    return reg_corner, ireg_corner, kappa

    