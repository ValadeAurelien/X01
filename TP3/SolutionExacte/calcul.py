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
       Matrix( [[2+sin(2*pi*x*epsilon), 0], [0, 4]] ),
       Matrix( [[2+sin(2*pi*x*epsilon), 0], [0, 4+sin(2*pi*x*epsilon)]] ),
       (2+sin(2*pi*x*epsilon))*(4+sin(2*pi*y*epsilon))* Matrix( [[1, 0], [0, 1]] ),
       (2+sin(10*x))* Matrix( [[1, 0], [0, 1]] )];

# --- Calc ---
for A in As:
    flux = A*gradu
    f = - (diff(flux[0], x) + diff(flux[1], y))
    print()
    print("*************************************************")
    print ("A :", A[0,0], A[0, 1], "\n   ",  A[1,0], A[1, 1])
    print()
    print("u :", u)
    print()
    print("gradu :", gradu[0], "\n       ", gradu[1])
    print()
    print("flux :", flux[0], "\n      ", flux[1])
    print()
    print("f :", str(simplify(f)).replace(")*c", ").*c").replace(")*s", ").*s").replace("**", "^").replace(" + (", " + ...\n    (").replace(" - (", " - ...\n    ("))
