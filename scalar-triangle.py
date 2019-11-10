import scipy as sc
from scipy.integrate import quad
from numpy import zeros,log,real
import matplotlib.pyplot as plt
import cmath
from  numpy.lib import scimath

def complex_quad(func,a1,a2):
    def Real(x):
        return sc.real(func(x))
    def Imag(x):
        return sc.imag(func(x))
    real_int = quad(Real, a1, a2)
    imag_int = quad(Imag, a1, a2)
    return (real_int[0],imag_int[0])

n = 10**(-6)
eps = -.00001
length = 300
m_pi = .100
m = .100 
s = sc.linspace(n,1,length)
Re = zeros((length))
Im = zeros((length))
t = .77526**2-1j*.77526*.1491


def a(y):
    return y
    
def b(x,y):
    return  y*(x-1.) 
    
def c(x,y,t):
    return (1.-x)*m**2 + x*(x-1.)*m**2 + x*t 
    
def d(x,y,t):
    return b(x,y)**2 - 4.*a(y)*c(x,y,t)
    
def ym(x,y,t):
    return (b(x,y)-scimath.sqrt(d(x,y,t)))/(2.*a(y))
def yp(x,y,t):
    return (b(x,y)+scimath.sqrt(d(x,y,t)))/(2.*a(y))


norm = complex_quad((lambda x: (1./scimath.sqrt(d(x,n,t))) * (scimath.log( (1.-x+ym(x,n,t))/ym(x,n,t) ) - scimath.log( (1.-x+yp(x,n,t))/yp(x,n,t) ))) ,0,1) 
print( norm)
    

for i in range(len(s)):
    e = s[i]
    val = complex_quad((lambda x: (1./scimath.sqrt(d(x,e,t))) * (scimath.log( (1.-x+ym(x,e,t))/ym(x,e,t) ) - scimath.log( (1.-x+yp(x,e,t))/yp(x,e,t) ))) ,0,1)
    Re[i] = val[0]
    Im[i] = val[1]

plt.xlabel(r'$\sqrt{s}/m_q$')
plt.ylabel(r'$\mathcal{M}(s)/\mathcal{M}(0)$')
plt.plot(scimath.sqrt(s)/m_pi,Re/norm[0],label='real part')
plt.plot(scimath.sqrt(s)/m,Im/norm[0],label='imaginary part')
plt.legend()
plt.show()
