import scipy as sc
from scipy.integrate import dblquad
from numpy import zeros,sqrt
import matplotlib.pyplot as plt
import cmath

def complex_quad(func,a1,b1,a2,b2):
    def Real(x,y):
        return sc.real(func(x,y))
    def Imag(x,y):
        return sc.imag(func(x,y))
    real_int = dblquad(Real,a1,b1,lambda x: a2,lambda x: b2)
    imag_int = dblquad(Imag,a1,b1,lambda x: a2,lambda x: b2)
    return (real_int[0],imag_int[0])

eps = 1
length = 15
m_pi = .135
m = .1
s = sc.linspace(4*m**2,0.1,length)
Re = zeros((length))
Im = zeros((length))
norm = complex_quad(lambda x1,x2: -1./(x1*x2*m**2 + x1*(1.-x1-x2)*m**2 - m**2 - 1j*eps),0.,1.,0.,1.)

for i in range(len(s)):
    a = complex_quad(lambda x1,x2: -1./(x1*x2*m**2 + x1*(1.-x1-x2)*m**2 + x2*(1.-x1-x2)*s[i] - m**2) ,0.,1.,0.,1.)
    #a = complex_quad(lambda x1,x2: -1./(x1*x2*m**2 + x1*(1.-x1-x2)*m**2 + x2*(1.-x1-x2)*s[i] - m**2 )+ sc.pi*1j*2/cmath.sqrt((1.-x1)**2 - 4.*(1.-x1+x1**2)*m**2/s[i]),0.,1.,0.,1.)
    #a = complex_quad(lambda x1,x2: -1./(x1*x2*m**2 + x1*(1.-x1-x2)*m**2 + x2*(1.-x1-x2)*s[i] - m**2 - 1j**eps),0.,1.,0.,1.)
    Re[i] = a[0]
    Im[i] = a[1]

plt.plot(sqrt(s)/m_pi,Re/norm[0],label='real part')
plt.plot(sqrt(s)/m_pi,Im/norm[1],label='imaginary part')
plt.legend()
plt.show()
