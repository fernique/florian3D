#####################################################
# Lemma 2 : R is a root of a quadratic polynomial P #
#####################################################

var('ab ac bc xc yc x y z R r_A r_B r_C')
assume(yc>0)
assume(z>0)
eqs=[
x^2+y^2+z^2-(R+r_A)^2,
(x-ab)^2+y^2+z^2-(R+r_B)^2,
(x-xc)^2+(y-yc)^2+z^2-(R+r_C)^2,
xc^2+yc^2-ac^2,
(xc-ab)^2+yc^2-bc^2
]

val={}

# get x
val.update({x:solve(eqs[0]==eqs[1],x)[0].right_hand_side()})
# get y
val.update({y:solve(eqs[0]==eqs[2],y)[0].right_hand_side().subs(val)})
# get xc
val.update({xc:solve(eqs[3]-eqs[4]==0,xc)[0].right_hand_side().subs(val)})
# get yc
val.update({yc:solve(eqs[4]==0,yc)[0].right_hand_side().subs(val)})


# subs in the first eq to get an eq P=0 in R and z (as well as ab,xc,yc,r_A,r_B,r_C)
P=eqs[0].subs(val).subs(val).numerator().expand()

# P is a quadratic polynomial A*R^2+B*R+C(z) with:
A=P.coefficient(R^2)
B=P.coefficient(R)
C=P.subs(R=0)

###################################################
# Lemma 3 : the coefficient C(z) of P is negative #
###################################################

# the derivative is easily seen to be negative:
factor(diff(C,z))

# now, let us show C(0)<0.
C0=(r_A^2+r_B^2-ab^2)*(r_A^2+r_C^2-ac^2)*(r_B^2+r_C^2-bc^2)
C1=-bc^4*r_A^2 + 2*bc^2*r_A^2*r_B^2 - r_A^2*r_B^4 + 2*bc^2*r_A^2*r_C^2 - r_A^2*r_C^4
C2=-ac^4*r_B^2 + 2*ac^2*r_B^2*r_C^2 - r_B^2*r_C^4 + 2*ac^2*r_B^2*r_A^2 - r_B^2*r_A^4
C3=-ab^4*r_C^2 + 2*ab^2*r_C^2*r_A^2 - r_C^2*r_A^4 + 2*ab^2*r_C^2*r_B^2 - r_C^2*r_B^4
C4=-2*r_A^2*r_B^2*r_C^2

# the following expression is 0
expand(C.subs(z=0)-C0-C1-C2-C3-C4)

# the following function is negative
var('x a b')
f(x)=x^4-(x^2-a^2)^2-(x^2-b^2)^2
factor(diff(f,x)(x))
expand(f(a+b))

###########################################################
# Lemma 4 : the smallest positive root of P is increasing #
###########################################################

# D is the discriminant of P
D=B^2-4*A*C

# D(0)>0 :
factor(D.subs(z=0))

# D=E*z^2+D(0)
# experiments show that E can be positive or not
# if E>0 then D>0 for any z, otherwise D>0 for z in some interval [z,z0]

# the numerator of the derivative of R_+ is positive
factor(diff((-B+sqrt(D))/2/A,z).numerator())

# the denominator is of the form sqrt(something complicate) hence positive also
diff((-B+sqrt(D))/2/A,z).denominator()
