module distaz_lib

  implicit none

  double precision, private, parameter :: pi = 3.14159265358979323846

contains
  subroutine distaz_cart(x1, y1, x2, y2, az, baz, dist)
    double precision, intent(in) :: x1, y1, x2, y2
    double precision, intent(out) :: az, baz, dist
    double precision :: dx, dy

    dx = x2 - x1
    dy = y2 - y1
    dist = sqrt(dx**2 + dy**2)

    az = atan2(dx, dy) * 180.0 / pi
    if (az < 0.0) then
      az = az + 360.0
    endif

    baz = az + 180.0
    if (baz >= 360.0) baz = baz - 360.0

  end subroutine distaz_cart

  subroutine distaz(stalat, stalon, evtlat, evtlon, az, baz, delta, dist)
    double precision, intent(in) :: stalat, stalon, evtlat, evtlon
    double precision, intent(out) :: delta, az, baz, dist
    double precision :: scolat, slon, ecolat, elon, deg2km
    double precision :: a,b,c,d,e,aa,bb,cc,dd,ee,g,gg,h,hh,k,kk
    double precision :: rhs1,rhs2,sph,rad,del,daz,dbaz,piby2, predel

    piby2=pi/2.
    rad=2.*pi/360.
    deg2km = 2*pi*6371/360

    sph=1.0/298.257
    scolat=piby2 - atan((1.-sph)*(1.-sph)*tan(stalat*rad))
    ecolat=piby2 - atan((1.-sph)*(1.-sph)*tan(evtlat*rad))
    slon=stalon*rad
    elon=evtlon*rad

    a=sin(scolat)*cos(slon)
    b=sin(scolat)*sin(slon)
    c=cos(scolat)
    d=sin(slon)
    e=-cos(slon)
    g=-c*e
    h=c*d
    k=-sin(scolat)

    aa=sin(ecolat)*cos(elon)
    bb=sin(ecolat)*sin(elon)
    cc=cos(ecolat)
    dd=sin(elon)
    ee=-cos(elon)
    gg=-cc*ee
    hh=cc*dd
    kk=-sin(ecolat)

    predel=a*aa + b*bb + c*cc
    if(abs(predel+1.).lt..0000000001) then
        predel=-1.
    endif
    if(abs(predel-1.).lt..0000000001) then
        predel=1.
    endif
    del=acos(predel)
    delta=del/rad
    dist=delta*deg2km

    rhs1=(aa-d)*(aa-d)+(bb-e)*(bb-e)+cc*cc - 2.
    rhs2=(aa-g)*(aa-g)+(bb-h)*(bb-h)+(cc-k)*(cc-k) - 2.
    dbaz=atan2(rhs1,rhs2)
    if(dbaz.lt.0.0d0) dbaz=dbaz+2*pi
    baz=dbaz/rad

    rhs1=(a-dd)*(a-dd)+(b-ee)*(b-ee)+c*c - 2.
    rhs2=(a-gg)*(a-gg)+(b-hh)*(b-hh)+(c-kk)*(c-kk) - 2.
    daz=atan2(rhs1,rhs2)
    if(daz.lt.0.0d0) daz=daz+2*pi
    az=daz/rad

    if(abs(baz-360.).lt..00001) baz=0.0
    if(abs(az-360.).lt..00001) az=0.0
  end subroutine distaz
end module distaz_lib
