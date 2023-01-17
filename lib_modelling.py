"""
module from AGII package
(c) Georg Kaufmann
"""

def grav_sphere(xo,yo,zo,xp,yp,zp,radius,rho):
    import numpy as np
    """
    !-----------------------------------------------------------------------
    ! Bouguer-Anomaly for a buried sphere
    ! x-axis is north, z-axis vertically down
    !
    ! Input parameters:
    ! xo,yo,zo [m]       - location of observation point
    ! xp,yp,zp [m]       - location of center of sphere
    ! rho [kg/m^3]       - density difference
    ! radius [m]         - radius of sphere
    !
    ! Output parameter:
    ! gx,gy,gz [mGal]    - gravity component at observation point
    !
    ! from:
    ! Blakeley (1995): Potential theory in gravity & magnetic applications
    !-----------------------------------------------------------------------
    """
    G       = 6.672e-11     # m^3/kg/s^2
    si2mgal = 1.e5          # m/s^2 -> mGal

    rx   = xo - xp
    ry   = yo - yp
    rz   = zo - zp
    r    = rx**2 + ry**2 + rz**2
    if (r == 0):
        print ('grav_sphere: obs. point in centrum of sphere')
    r    = np.sqrt(r)
    r3   = r**3
    mass = 4./3.*np.pi*rho*radius**3
    gx = -G * mass * rx/ r3 * si2mgal
    gy = -G * mass * ry/ r3 * si2mgal
    gz = -G * mass * rz/ r3 * si2mgal
    return gx,gy,gz


def mag_sphere (xo,yo,zo,xp,yp,zp,radius,earthincl,earthdecl,mag,incl,decl,theta):
    import numpy as np
    """
    !-----------------------------------------------------------------------
    ! Magnetic anomaly of a uniformly magnetised sphere.
    ! x-axis is north, z-axis vertically down
    !
    ! Input parameters:
    ! xo,yo,zo [m]         - location of observation point
    ! xp,yp,zp [m]         - location of center of sphere
    ! radius [m]           - radius of sphere
    ! earthincl [degrees]  - inclination of earth field, positive below horizontal
    ! earthdecl [degrees]  - declination of earth field, positive east of true north
    ! mag [A/m]            - magnetisation      
    ! incl [degrees]       - inclination of magnetisation, positive below horizontal
    ! decl [degrees]       - declination of magnetisation, positive east of true north
    ! theta [degrees]      - azimuth of x axis in degrees positive east to north
    !
    ! Output parameter:
    ! bx,by,bz [nT]        - Magnetic induction component at observation point
    ! t [nT]               - Total field at observation point
    !
    ! from:
    ! Blakeley (1995): Potential theory in gravity & magnetic applications
    !-----------------------------------------------------------------------
    """
    mu04pi  = 1.000e-7      # Vs/Am
    t2nt    = 1.e9          # T -> nT
    mx,my,mz = dircos (incl,decl,theta)
    fx,fy,fz = dircos (earthincl,earthdecl,theta)

    rx   = xo - xp
    ry   = yo - yp
    rz   = zo - zp
    r2   = rx**2 + ry**2 + rz**2
    if (r2 == 0):
        print ('mag_sphere: obs point in centrum of sphere')
    r    = np.sqrt(r2)
    r5   = r**5
    dot  = rx*mx + ry*my + rz*mz
    moment = 4./3.*np.pi*mag*radius**3
    bx = t2nt * mu04pi * moment * (3.*dot*rx - r2*mx) / r5
    by = t2nt * mu04pi * moment * (3.*dot*ry - r2*my) / r5
    bz = t2nt * mu04pi * moment * (3.*dot*rz - r2*mz) / r5
    t  = fx*bx + fy*by + fz*bz
    return bx,by,bz,t


def dircos(incl,decl,azim):
    import numpy as np
    """
    !-----------------------------------------------------------------------
    ! Input parameters:
    ! incl [degrees]     - inclination positive below horizontal
    ! decl [degrees]     - declination positive east from north
    ! azim [degrees]     - azimuth of x axis in degrees positive east to north
    ! Output parameter:
    ! a,b,c              - three direction cosines
    !
    ! from: 
    ! Blakeley (1995): Potential theory in gravity & magnetic applications
    !-----------------------------------------------------------------------
    """
    deg2rad = np.pi / 180.

    a = np.cos(incl*deg2rad) * np.cos((decl-azim)*deg2rad)
    b = np.cos(incl*deg2rad) * np.sin((decl-azim)*deg2rad)
    c = np.sin(incl*deg2rad)
    return a,b,c
