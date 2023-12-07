K.<x,y,mu> = QQ[]
#In our mind: -1 \leq y \leq 0 ; y\leq x \leq -y.

# Basics
e1 = vector(QQ, [1,0])
e2 = vector(QQ, [0,1])
eins = e1+e2

#Vertices of a reduced triangle T
v0 = -eins
v1 = x*e1 + e2
v2 = e1 + y*e2

#Let Ei be the edge of T opposite to vi
#Computing the equations for Ei
m0 = Matrix([v1,v2])
a0 = m0.solve_right(eins)

#x in aff(E0) iff a0*x = 1
#aff(E0) = {(u,v) : (y-1)*u + (x-1)*v = x*y-1}

m1 = Matrix([v0,v2])
a1 = m1.solve_right(eins)

#aff(E1) = {(u,v) : (y+1)*u -2*v = 1-y}

m2 = Matrix([v0,v1])
a2 = m2.solve_right(eins)

#aff(E2) = {(u,v) : -2*u + (x+1)*v = 1-x}

#point on E2 to be varied using mu
p2 = (1-mu)*v0 + mu*v1

#the other two points that could form a rectangle with p2 so that p2 is the corner
A0 = Matrix([a0,e2])
p0 = A0.solve_right(e1 + p2[1]*e2)

A1 = Matrix([a1, e1])
p1 = A1.solve_right(e1 + p2[0]*e2)

#squared edge lengths of the two catheters
d0 = (p0-p2).dot_product(p0-p2)
d1 = (p1-p2).dot_product(p1-p2)

#Finding the mu for which d0=d1
objective = (d0-d1)
objective = objective.numerator()

c2 = objective.coefficient(mu^2)
c1 = objective.coefficient(mu)
c0 = objective(mu=0)

p=c1/c2
q=c0/c2

sol1 = -p/2 + sqrt((p/2)^2 -q)
sol2 = -p/2 - sqrt((p/2)^2 -q)

#sol1 is always valid, sol2 never

sol111 = -2/(y-3) #this is actually sol1, but nicer....

hurkens_function = d0(mu=sol111)
hurkens_gradient1 = hurkens_function.derivative(x).factor()
hurkens_gradient2 = hurkens_function.derivative(y).factor()

#One sees that there are no critical points in the interior of the domain.

#It remains to check the boundary.

bfunction0 = hurkens_function(y=-1) #identically 1
bfunction_minus = hurkens_function(x=y)
bfunction_plus = hurkens_function(x=-y)

#The following is the unique critical point of both bfunction_plus and .._minus in (-1,0).

critical_y = QQbar(3-2*sqrt(3))

#Let's evaluate...

cand_m = bfunction_minus(y=critical_y)
cand_p = bfunction_plus(y=critical_y)

#We see that cand_p < cand_m. So cand_p realizes the minimum.

#Recall that cand_p = hurkens_function(-critical_y,critical_y).
#Thus, 

flt2 = 2/sqrt(cand_p)


minimal_poly = flt2.minpoly()

# This shows that flt2 = 1+2/sqrt(3).
