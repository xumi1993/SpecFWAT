subroutine distaz(stalat, stalon, evtlat, evtlon, az, baz, delta, dist)
    implicit none
    double precision :: stalat, stalon, evtlat, evtlon, delta, az, baz, dist
    double precision :: scolat, slon, ecolat, elon, deg2km
    double precision :: a,b,c,d,e,aa,bb,cc,dd,ee,g,gg,h,hh,k,kk
    double precision :: rhs1,rhs2,sph,rad,del,daz,dbaz,pi,piby2, predel

    pi=3.141592654
    piby2=pi/2.
    rad=2.*pi/360.
    deg2km = 2*pi*6371/360

    sph=1.0/298.257
    scolat=piby2 - atan((1.-sph)*(1.-sph)*tan(dble(stalat)*rad))
    ecolat=piby2 - atan((1.-sph)*(1.-sph)*tan(dble(evtlat)*rad))
    slon=dble(stalon)*rad
    elon=dble(evtlon)*rad

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
