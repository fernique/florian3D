# the edge lengths of the tetrahedron, the radii of its spheres and the radius R of its support sphere satisfies a polynomial relation P=0 that we compute here
# good point: P turns out to be quadratic in R (as well as in the square of each edge length)
# bad point: P is the sum of 420 monomials each of degree 8

# edge lengths
var('ab ac ad bc bd cd')
# radii (R is the radius of the support sphere)
var('ra rb rc rd R')

# coordinates of the vertices A,B,C,D, as well as (x,y,z) the center of the support sphere
var('xa ya za xb yb zb xc yc zc xd yd zd x y z')

# all these quantities must be equal 0 -> system of polynomial equations
eqs=[
(xb-xa)^2+(yb-ya)^2+(zb-za)^2-ab^2,
(xc-xa)^2+(yc-ya)^2+(zc-za)^2-ac^2,
(xd-xa)^2+(yd-ya)^2+(zd-za)^2-ad^2,
(xc-xb)^2+(yc-yb)^2+(zc-zb)^2-bc^2,
(xd-xb)^2+(yd-yb)^2+(zd-zb)^2-bd^2,
(xd-xc)^2+(yd-yc)^2+(zd-zc)^2-cd^2,
(xa-x)^2+(ya-y)^2+(za-z)^2-(ra+R)^2,
(xb-x)^2+(yb-y)^2+(zb-z)^2-(rb+R)^2,
(xc-x)^2+(yc-y)^2+(zc-z)^2-(rc+R)^2,
(xd-x)^2+(yd-y)^2+(zd-z)^2-(rd+R)^2
]

# dictionnary of known values
val={}

# we set the axis so that A=(0,0,0), B=(ab,0,0), C=(xc,yc,0), yc>0 and zd>0
val.update({xa:0,ya:0,za:0,xb:ab,yb:0,zb:0,zc:0})
assume(yc>0)
assume(zd>0)

# update equations
eqs=[i.subs(val) for i in eqs]

# first get xc
val.update({xc:solve(eqs[1]==eqs[3],xc)[0].right_hand_side()})
# then yc
val.update({yc:solve(eqs[1]==0,yc)[0].right_hand_side().subs(val)})
# xd
val.update({xd:solve(eqs[2]==eqs[4],xd)[0].right_hand_side()})
# then yd
val.update({yd:solve(eqs[2]==eqs[5],yd)[0].right_hand_side().subs(val)})
# then zd
val.update({zd:solve(eqs[2]==0,zd)[0].right_hand_side().subs(val)})
# now, use the last 4 equations to get a linear system in x,y,z,R (and the other already determined variables)
# solve it to get x,y,z
t=solve([eqs[6]==eqs[7],eqs[6]==eqs[8],eqs[6]==eqs[9]],[x,y,z])[0]
val.update({
x:t[0].right_hand_side().subs(val).full_simplify(),
y:t[1].right_hand_side().subs(val).full_simplify(),
z:t[2].right_hand_side().subs(val).full_simplify()
})

# injecting in eqs[6] which has variables x,y,z,R yields the wanted polynomial P
P=eqs[6].subs(val).numerator().expand()


