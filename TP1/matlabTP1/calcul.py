##
from sympy import *

x = Symbol('x')
y = Symbol('y')
b = Symbol('alpha')

u = cos(pi*x)*cos(2*pi*y)
A = sin(b*pi*x)*sin(b*pi*y)+2
#A = 1

f = u - diff(A, x) * diff(u, x) + diff(A, y) * diff(u, y) - A *( diff(diff(u, x), x) + diff(diff(u, y), y))

#print(f)
print(str(simplify(f)).replace(")*c", ").*c").replace(")*s", ").*s").replace("**", "^"))

##

from sympy import *
#A = [ [ a , b ], \
#      [ c , d ]]
#telle que A(X-B) renvoie sur le triangle unité 

a = Symbol('a')
b = Symbol('b')
c = Symbol('c')
d = Symbol('d')

x1 = Symbol('x1')
y1 = Symbol('y1')
x2 = Symbol('x2')
y2 = Symbol('y2')
x3 = Symbol('x3')
y3 = Symbol('y3')

ToSolve = Matrix( ( (x2-x1, y2-y1, 0, 0, 1), (0, 0, x2-x1, y2-y1, 0), (x3-x1, y3-y1, 0, 0, 0), (0, 0, x3-x1, y3-y1, 1) ) )
system = ToSolve[:,:-1], ToSolve[:,-1]
A = linsolve(system, a, b, c, d)
print(XX)

