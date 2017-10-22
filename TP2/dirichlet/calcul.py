from sympy import *

x = Symbol('x')
y = Symbol('y')
epsilon = Symbol('epsilon')

# --- U ---
u = sin(pi*x)*sin(2*pi*y)
gradu = Matrix( [[diff(u, x)], [diff(u, y)]] )

# --- A ---
As = [ Matrix( [[1, 0], [0, 1]]),
       (sin(alpha*pi*x)*sin(alpha*pi*y)+2) * Matrix( [[1, 0], [0, 1]]) ]

# --- Calc ---
for A in As:
    flux = A*gradu
    f = u - ( diff(flux[0], x) + diff(flux[1], y) )
    print()
    print (A)
    print(str(simplify(f)).replace(")*c", ").*c").replace(")*s", ").*s").replace("**", "^"))


