

# P is the polynomial relation between edge length, radii and radius of the support sphere
# it is quadratic in R^2 (as well as in the square of each edge length)
# it is defined in the file radius.sage (as well as the variables ab,ac,ad,bc,bd,cd,R)


#####################################################################
# The point is to simplify the formulas for the coefficients        #
# Depending on the cases after dimension reduction                  #
# We heuristically minimize the number of occurrences of variables  #
#####################################################################

var('r')

# count the nb of occurrences of variables in the SR exp
score = lambda exp: sum([str(exp).count(str(x)) for x in exp.variables() if str(x)!='r'])

# for the sake of lisibility we use xi instead of square of edge length
var('x0 x1 x2 x3 x4 x5')
rename={
ab^2:x0,ac^2:x1,ad^2:x2,bc^2:x3,bd^2:x4,cd^2:x5,
ab^4:x0^2,ac^4:x1^2,ad^4:x2^2,bc^4:x3^2,bd^4:x4^2,cd^4:x5^2,
}

Z=NumberField((r+1)^2-2,'r')
#Z.inject_variables()
K=PolynomialRing(Z,'x',6)
#K.inject_variables()


# heuristic to find an expression equivalent to a given SR(polynomial) with as less variable occurences as possible
def simplify(e):
    if e.degree()<2:
        return str(e)

    # TODO : optimize further when e.degree()==2 ?
    # a*x*y+b*x+c*y -> (x+c/a)*(a*y+b)-b*c/a
    # a*x^2+b*x*y+c*y^2 -> a*(x+b/2/a*y)^2+(c-b^2/4/a)*y^2      (for y=1 as well)

    # if factorizable then apply recursively on each factor
    if not e.is_prime():
        c=e.factor().unit()
        s='' if c==1 else '-' if c==-1 else str(c)+'*'
        for (f,d) in e.factor():
            f=simplify(f)
            s+='('+f+')'
            if d>1:
                s=s+'^'+str(d)
            s+='*'
        s=s[0:-1] ## remove the last *
        return s

    # otherwise, find x^d which maximizes the nb of occurences in e
    L=[]
    for x in e.variables():
        for d in range(1,e.degree(x)+1):
            L.append((e.coefficient(x^d).number_of_terms(),x^d))
    L.sort(reverse=True)
    xd=L[0][1]
    
    # then factor e=f*x^d+g
    (f,g)=(0,0)
    for (c,m) in list(e):
        if (xd).divides(m):
            f+=K(c*m/xd)
        else:
            g+=c*m

    # and apply recursively
    g=simplify(g)
    f=str(xd)+'*('+simplify(f)+')'
    if g[0]!='-':
        f+='+'
    return f+g



#####################
# conversion to cpp #
#####################

# convert SR(polynomial) to a C++ string
def cpp(exp):
    (x,y)=(SR.wild(0),SR.wild(1))
    if exp.is_integer():
        return 'I('+str(exp)+')' if exp>0 else '-I('+str(-exp)+')'
    if exp in QQ:
        return 'I('+str(exp.numerator())+')/I('+str(exp.denominator())+')'
    if exp.is_symbol():
        return str(exp)
    m=exp.match(x*y)
    if m !=None:
        (X,Y)=(m[y],m[x]) if m[y].is_integer() else (m[x],m[y])
        if X==1:
            return cpp(Y)
        if X==-1:
            return '-'+cpp(Y)
        X=cpp(X) if X.is_integer() or X.is_symbol() else '('+cpp(X)+')'
        Y=cpp(Y) if Y.is_integer() or Y.is_symbol() else '('+cpp(Y)+')'
        return X+'*'+Y
    m=exp.match(x/y)
    if m !=None:
        (X,Y)=(m[y],m[x]) if m[y].is_integer() else (m[x],m[y])
        if X==1:
            return cpp(Y)
        if X==-1:
            return '-'+cpp(Y)
        X=cpp(X) if X.is_integer() or X.is_symbol() else '('+cpp(X)+')'
        Y=cpp(Y) if Y.is_integer() or Y.is_symbol() else '('+cpp(Y)+')'
        return X+'/'+Y
    m=exp.match(x+y)
    if m !=None:
        X=cpp(m[x])
        Y=cpp(m[y])
        if X[0]=='-' and Y[0]!='-':
            return Y+X
        if X[0]!='-' and Y[0]=='-':
            return X+Y
        if X[0]=='-' and Y[0]=='-':
            return X+Y
        if X[0]!='-' and Y[0]!='-':
            return X+'+'+Y
    m=exp.match(x^2)
    if m !=None:
        return 'square('+cpp(m[x])+')'
    m=exp.match(sqrt(x))
    if m !=None:
        return 'sqrt('+cpp(m[x])+')'
    m=exp.match(x^y)
    if m !=None:
        return 'pow('+cpp(m[x])+','+str(m[y])+')'
    print(exp)
    raise ValueError # mismatch (should not happen)

#####################
# let's go optimize #
#####################

# first define A,B,C,D to be optimized

# to compute the radius R from the edge lengths, we see P as a quadratic polynomial in R
A=P.coefficient(R^2).expand()
B=P.coefficient(R).expand()
C=P.subs(R=0).expand()
D=B^2-4*A*C
# remark: the cayley_menger determinant (288*V^2) divides D


# now we optimize


# dictionnary to replace radius and edges along which contact is assumed (or radius of the support sphere) by their value
case1={ra:1,rb:1,rc:1,rd:1,ab:2}
case2={ra:1,rb:1,rc:1,rd:r,ab:2}
case3={ra:1,rb:1,rc:1,rd:r,ad:1+r}
case4={ra:1,rb:1,rc:r,rd:r,ab:2}
case5={ra:1,rb:1,rc:r,rd:r,ad:1+r}
case6={ra:1,rb:1,rc:r,rd:r,cd:r+r}
case7={ra:1,rb:r,rc:r,rd:r,ad:1+r}
case8={ra:1,rb:r,rc:r,rd:r,cd:r+r}
case9={ra:r,rb:r,rc:r,rd:r,cd:r+r}

(s1,s2)=(0,0)

dico=case1

[score(i) for i in [A,B,C,D]]

A2=A.subs(dico)
B2=B.subs(dico)
C2=C.subs(dico)
D2=D.subs(dico)

A2=expand(A2.subs(rename))
B2=expand(B2.subs(rename))
C2=expand(C2.subs(rename))
D2=expand(D2.subs(rename))

[score(i) for i in [A2,B2,C2,D2]]
# case6 D : 2661 var ! -> 45
s1+=sum([score(i) for i in [A2,B2,C2,D2]])

A2=SR(simplify(K(str(A2))))
B2=SR(simplify(K(str(B2))))
C2=SR(simplify(K(str(C2))))
D2=SR(simplify(K(str(D2))))

[score(i) for i in [A2,B2,C2,D2]]
s2+=sum([score(i) for i in [A2,B2,C2,D2]])

print('\tI A='+cpp(A2)+';\n'+'\tI B='+cpp(B2)+';\n'+'\tI C='+cpp(C2)+';\n'+'\tI D='+cpp(D2)+';\n')

