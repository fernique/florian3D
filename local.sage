#'SageMath version 8.6, Release Date: 2019-01-15'

# We here compute, for each tetrahedron claimed to be the densest in Theorem 1, an explicit neighborhood over which this holds


# we assume P is the polynomial defined in radius_and_ac.sage
var('r')
V=var('ab ac ad bc bd cd')
var('ra rb rc rd R')

############
# routines #
############

# cosine of angle opposite to a in triangle with edge length a, b, and c (cosine law)
cosa=lambda a,b,c: (b/c+c/b-a^2/b/c)/2

# cosine of the dihedral angle along AB (spherical law of cosines)
#cosdi=lambda ab,ac,ad,bc,bd,cd: (cosa(cd,ac,ad)-cosa(bd,ab,ad)*cosa(bc,ab,ac))/sqrt(((1-cosa(bd,ab,ad)^2)*(1-cosa(bc,ab,ac)^2)))

# volume of the tetrahedron
#vol=lambda ab,ac,ad,bc,bd,cd: 1/6*ab*ac*ad*sqrt((1+2*cosa(bd,ab,ad)*cosa(cd,ac,ad)*cosa(bc,ab,ac)-cosa(bd,ab,ad)^2-cosa(cd,ac,ad)^2-cosa(bc,ab,ac)^2))
# full_simplify of the previous expression (better precision in RIF)
#vol=lambda ab,ac,ad,bc,bd,cd: 1/12*sqrt(-(ad^2*bc^4 + ac^2*bd^4 + ab^2*cd^4 + (ab^2*ac^2 + ad^4 - (ab^2 + ac^2)*ad^2)*bc^2 - (ab^2*ac^2 - ac^4 - (ab^2 - ac^2)*ad^2 + (ac^2 + ad^2)*bc^2)*bd^2 + (ab^4 - ab^2*ac^2 - (ab^2 - ac^2)*ad^2 - (ab^2 + ad^2)*bc^2 - (ab^2 + ac^2 - bc^2)*bd^2)*cd^2))
# another equivalent formula homemade simplification
vol=lambda ab,ac,ad,bc,bd,cd: 1/12*sqrt(-((ac + ad)*(ac - ad)*(bc + bd)*(bc - bd) + (ab^2 - ac^2 - ad^2 - bc^2 - (bd + cd)*(bd - cd))*cd^2)*ab^2 + (ad^2*(bc + cd)*(bc - cd) - (ac^2 - ad^2 - bc^2 + (bd + cd)*(bd - cd))*bd^2)*ac^2 - (bd^2*cd^2 + (ad^2 + bc^2 - bd^2 - cd^2)*ad^2)*bc^2)
# ATTENTION : the eventual value of epsilon is highly dependent on used the formula !!! we used the last one

# solid angle in A (Lagrange's theorem)
def solid(ab,ac,ad,bc,bd,cd):
    # f: denominator in the Lagrange formula
    f=(2*ab*ac*ad+(ac^2+ad^2-cd^2)*ab+(ab^2+ad^2-bd^2)*ac+(ab^2+ac^2-bc^2)*ad)/12
    V=vol(ab,ac,ad,bc,bd,cd)
    A=2*arctan(V/f) # solid angle in [-pi,pi]
    if parent(A)==RIF and A.upper()<=0:
        A+=2*pi # solid angle in [0,2pi]
    return A

# permutations on the edge lengths yielded by an isometry which moves A,B,C,D (resp.) on A
permut=lambda ab,ac,ad,bc,bd,cd: [(ab,ac,ad,bc,bd,cd),(bc,bd,ab,cd,ac,ad),(cd,ac,bc,ad,bd,ab),(ad,bd,cd,ab,ac,bc)]

# volume of the part covered by balls
cov=lambda e,rad: vector(solid(*i) for i in permut(*e))*vector([i^3/3 for i in rad])

# density
density=lambda edges,radii: cov(edges,radii)/vol(*edges)

###############
# T1111 tight #
###############

radii={ra:1,rb:1,rc:1,rd:1}

# eps-neighborhood (the edges can only be elongated)
eps=1/46
edges={
ab:ra+rb+RIF(0,eps),
ac:ra+rc+RIF(0,eps),
ad:ra+rd+RIF(0,eps),
bc:rb+rc+RIF(0,eps),
bd:rb+rd+RIF(0,eps),
cd:rc+rd+RIF(0,eps)
}

# the derivative of the density are all negatives over this neighborhood
# by symmetry, one suffices
[RIF(diff(density((ab,ac,ad,bc,bd,cd),(ra,rb,rc,rd)),v).subs(edges).subs(radii)).upper() for v in V]

# Remarks:
# 1) the derivatives yield different intervals...because the symmetry of the formula is broken by optimization
# 2) the value of epsilon is highly dependent on the expression of the derivatives...(which may depends on the sage version)

###############
# T11rr tight #
###############

radii={ra:1,rb:1,rc:r,rd:r}

# eps-neighborhood (the edges can only be elongated)
eps=1/203
edges={
ab:ra+rb+RIF(0,eps),
ac:ra+rc+RIF(0,eps),
ad:ra+rd+RIF(0,eps),
bc:rb+rc+RIF(0,eps),
bd:rb+rd+RIF(0,eps),
cd:rc+rd+RIF(0,eps)}

# the derivative of the density are all negatives over this neighborhood
# by symmetry A/B and C/D: would be sufficient to check for ab, cd and one of the 4 other ones
[RIF(diff(density((ab,ac,ad,bc,bd,cd),(ra,rb,rc,rd)),v).subs(edges).subs(radii).subs(r=RIF(sqrt(2)-1))).upper() for v in V]


###############
# T1rrr tight #
###############

radii={ra:1,rb:r,rc:r,rd:r}

# eps-neighborhood (the edges can only be elongated)
eps=1/148
edges={
ab:ra+rb+RIF(0,eps),
ac:ra+rc+RIF(0,eps),
ad:ra+rd+RIF(0,eps),
bc:rb+rc+RIF(0,eps),
bd:rb+rd+RIF(0,eps),
cd:rc+rd+RIF(0,eps)}

# the derivative of the density are all negatives over this neighborhood
[RIF(diff(density((ab,ac,ad,bc,bd,cd),(ra,rb,rc,rd)),v).subs(edges).subs(radii).subs(r=RIF(sqrt(2)-1))).upper() for v in V]




######################
# Trrrr bi-stretched #
######################

radii={ra:r,rb:r,rc:r,rd:r}

# we first compute formulas of the partial derivatives of the radius of the support sphere
# P is the polynomial defined in radius.sage
A=P.coefficient(R^2).expand()
B=P.coefficient(R).expand()
C=P.subs(R=0).expand()
R0=(-B+sqrt(B^2-4*A*C))/2/A # solve(P==0,R)[0].right_hand_side() (worse formula in practice)
dR0=vector([diff(R0,v) for v in V])

# define the function f s.t. f(T^*)=0 
s=r*sqrt(2*sqrt(6) + 6) # length of the two stretched edges rr in Trrrr^*, say bd and cd
d0=density((ra+rb,ra+rc,ra+rd,rb+rc,s,s),(ra,rb,rc,rd))
a=r/d0 # yields eps=1/2963
a=r/d0*((ra+rb)/ab)^35*((ra+rc)/ac)^35*((ra+rd)/ad)^58*((rb+rc)/bc)^40*(s/bd)^35*(s/cd)^35 # 1/173
# heuristic to fix the exponents: fix some "bad" eps and iteratively increase the exponent which corresponds to the most negative derivative (while they are not all positive) ; or even choose the exponent proportional to the negativity of the derivative
# we shall not use f but only df
# f=(R0(ab,ac,ad,bc,bd,cd)-a*density((ab,ac,ad,bc,bd,cd),(ra,rb,rc,rd)))

# partial derivative must be nonnegative over the neighborhood
df=dR0-vector([diff(a*density((ab,ac,ad,bc,bd,cd),(ra,rb,rc,rd)),v) for v in V])

# search for a neighborhood (both stretched edges bd and cd can be elongated or reduced)
edges=lambda eps:{
ab:ra+rb+RIF(0,eps),
ac:ra+rc+RIF(0,eps),
ad:ra+rd+RIF(0,eps),
bc:rb+rc+RIF(0,eps),
bd:s+RIF(-eps,eps),
cd:s+RIF(-eps,eps)
}

# these quantities must all be positive
[RIF(i).lower() for i in df.subs(edges(1/173)).subs(radii).subs(r=RIF(sqrt(2)-1))]


######################
# T111r 11-stretched #
######################

radii={ra:r,rb:1,rc:1,rd:1}

# we first compute formulas of the partial derivatives of the radius of the support sphere
A=P.coefficient(R^2).expand()
B=P.coefficient(R).expand()
C=P.subs(R=0).expand()
R0=(-B+sqrt(B^2-4*A*C))/2/A # solve(P==0,R)[0].right_hand_side() (worse formula in practice)
dR0=vector([diff(R0,v) for v in V])


# define the function f s.t. f(T^*)=0 
t=4*sqrt(2)*sqrt(r/(2*r+1)) # length of the stretched 11-edge, say CD
d0=density((ra+rb,ra+rc,ra+rd,rb+rc,rb+rd,t),(ra,rb,rc,rd))
a=r/d0 # fails
a=r/d0*((ra+rb)/ab)^2*((ra+rc)/ac)^8*((ra+rd)/ad)^8*((rb+rc)/bc)^4*((rb+rd)/bd)^4*(t/cd)^16 # 1/2008
a=r/d0*((ra+rb)/ab)^20*((ra+rc)/ac)^80*((ra+rd)/ad)^80*((rb+rc)/bc)^40*((rb+rd)/bd)^40*(t/cd)^160 # 1/445
# a=r/d0*((ra+rb)/ab)^200*((ra+rc)/ac)^800*((ra+rd)/ad)^800*((rb+rc)/bc)^400*((rb+rd)/bd)^400*(t/cd)^1600 # 1/714
# we shall not use f but only df
# f=(R0(ab,ac,ad,bc,bd,cd)-a*density((ab,ac,ad,bc,bd,cd),(ra,rb,rc,rd)))

# partial derivative must be nonnegative over the neighborhood
df=dR0-vector([diff(a*density((ab,ac,ad,bc,bd,cd),(ra,rb,rc,rd)),v) for v in V])

# search for a neighborhood (the stretched 11-edge ac can be elongated or reduced)
edges=lambda eps:{
ab:ra+rb+RIF(0,eps),
ac:ra+rc+RIF(0,eps),
ad:ra+rd+RIF(0,eps),
bc:rb+rc+RIF(0,eps),
bd:rb+rd+RIF(0,eps),
cd:t+RIF(-eps,eps)
}

# these quantities must all be positive
[RIF(i).lower() for i in df.subs(edges(1/445)).subs(radii).subs(r=RIF(sqrt(2)-1))]
