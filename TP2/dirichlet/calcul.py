##
from sympy import *

x = Symbol('x')
y = Symbol('y')
b = Symbol('alpha')

u = sin(pi*x)*sin(2*pi*y)
A = sin(b*pi*x)*sin(b*pi*y)+2
#A = 1

f = u - diff(A, x) * diff(u, x) + diff(A, y) * diff(u, y) - A *( diff(diff(u, x), x) + diff(diff(u, y), y))

#print(f)
print(str(simplify(f)).replace(")*c", ").*c").replace(")*s", ").*s").replace("**", "^"))


