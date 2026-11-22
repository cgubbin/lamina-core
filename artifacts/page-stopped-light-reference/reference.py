import cmath
def residual(w,b):
    k=w/299.792458
    metal=4-3130**2/(w*(w+10.7j))
    def root(e):
        s=cmath.sqrt(e*k*k-b*b)
        return s if s.imag>0 else -s
    h,d=1,-1j*root(metal)/metal
    for e,t in [(11.68,.110),(metal,.050)]:
        s=root(e)
        c,a=cmath.cos(s*t),cmath.sin(s*t)
        h,d=c*h+e*a/s*d,-s*a/e*h+c*d
    return d-1j*root(1)*h
def solve(b,w):
    for i in range(30):
        f=residual(w,b)
        dw=f/((residual(w+.001,b)-residual(w-.001,b))/.002)
        w-=dw
        if abs(dw)<1e-10:return w
    raise ArithmeticError
for b,w in [(7.41,805-4.14j),(30.3,778-5.06j)]:
    print(b,solve(b,w),(solve(b+.001,w)-solve(b-.001,w))/.002)
# Independent validation of the passive TM branch. Distances are in micrometres.
from pathlib import Path
import csv
out=Path(__file__).parent
with (out/'dispersion.csv').open('w') as f:
    writer=csv.writer(f)
    writer.writerow(['beta_per_um','omega_real_per_ps','omega_imag_per_ps'])
    w=805-4j
    for i in range(331):
        b=4+i/10
        w=solve(b,w)
        writer.writerow([b,w.real,w.imag])
with (out/'fields.csv').open('w') as f:
    writer=csv.writer(f)
    writer.writerow(['beta_per_um','depth_um','abs_Hy','abs_Ex_scaled','abs_Ez_scaled'])
    for b in (6,7.41,9,26,30.3,34):
        w=solve(b,805-4j)
        k=w/299.792458
        eps=4-3130**2/(w*(w+10.7j))
        def transverse(e):
            s=cmath.sqrt(e*k*k-b*b)
            return s if s.imag>0 else -s
        for i in range(-100,301):
            z=i/1000
            h,d=1,-1j*transverse(eps)/eps
            e=eps
            if z<0:
                factor=cmath.exp(-1j*transverse(eps)*z)
                h,d=h*factor,d*factor
            else:
                for e,start,end in [(11.68,0,.11),(eps,.11,.16),(1,.16,1)]:
                    s=transverse(e)
                    t=min(z,end)-start
                    c,a=cmath.cos(s*t),cmath.sin(s*t)
                    h,d=c*h+e*a/s*d,-s*a/e*h+c*d
                    if z<=end:break
            writer.writerow([b,z,abs(h),abs(d/k),abs(b*h/(k*e))])

# Locate extrema independently of the thesis's rounded target wavevectors.
def real_slope(b, w):
    step=0.0001
    return ((solve(b+step,w)-solve(b-step,w))/(2*step)).real

for lo,hi,w in [(7,8,805-4j),(29,32,778-5j)]:
    sign=real_slope(lo,w)
    assert sign*real_slope(hi,w)<0
    for _ in range(35):
        mid=(lo+hi)/2
        if real_slope(mid,w)*sign>0:
            lo=mid
        else:
            hi=mid
    print("stationary beta and omega:",mid,solve(mid,w))
