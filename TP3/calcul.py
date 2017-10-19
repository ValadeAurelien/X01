from sympy import *

x = Symbol('x')
y = Symbol('y')
epsilon = Symbol('epsilon')

# --- U ---
u = sin(pi*x)*sin(pi*y)
gradu = Matrix( [[diff(u, x)], [diff(u, y)]] )

# --- A ---
As = [ Matrix( [[1, 0], [0, 1]]),
       Matrix( [[1, 0], [0, 2]]),
       Matrix( [[2+2*sin(2*pi*x/epsilon), 0], [0, 4]] ),
       Matrix( [[2+2*sin(2*pi*x/epsilon), 0], [0, 4+sin(2*pi*x/epsilon)]] ),
       (2+2*sin(2*pi*x/epsilon))*(4+sin(2*pi*x/epsilon))* Matrix( [[1, 0], [0, 1]] ) ]

# --- Calc ---
for A in As:
    flux = A*gradu
    f = diff(flux[0], x) + diff(flux[1], y)
    print()
    print (A)
    print(str(simplify(f)).replace(")*c", ").*c").replace(")*s", ").*s").replace("**", "^"))


