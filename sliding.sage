
# claimed maximal density (approximation from below)
d1111=0.779635570044251
d111r=0.812542027810809
d11rr=0.810466032832061
d1rrr=0.806503318194779
drrrr=0.784688454045208

# used variables
var('ra rb rc d r')

# the lower bound on the volume of P (divided by h)
# i.e. the left hand side in the inequality of the sliding Lemma
x=(rb+rc)*(ra+(1-d)*max(rb,rc))^2/6/(1+r)

# the upper bound on the volume of S_B\cap S_C (divided by h and then for the limit h=0)
# i.e. the right hand side in the inequality of the sliding Lemma, using the minoration of 
y=d*min(rb,rc)*(d-1)^2*(rb+rc)^2/2/min(ra,rb,rc)

# we have to check that z:=x-y is positive for each face XYZ of a tetrahedron of type t (the points towards which we slide plays no role)
z=x-y

# tetrahedron of type 1111: face 111
RIF(z.subs({ra:1,rb:1,rc:1,d:d1111^(1/3)}).subs(r=sqrt(2)-1)).lower() # sliding 1 in face 111

# tetrahedron of type 111r: faces 111, 11r
RIF(z.subs({ra:1,rb:1,rc:1,d:d111r^(1/3)}).subs(r=sqrt(2)-1)).lower() # sliding 1 in face 111
RIF(z.subs({ra:1,rb:1,rc:r,d:d111r^(1/3)}).subs(r=sqrt(2)-1)).lower() # sliding 1 in face 11r
RIF(z.subs({ra:r,rb:1,rc:1,d:d111r^(1/3)}).subs(r=sqrt(2)-1)).lower() # sliding r in face 11r

# tetrahedron of type 11rr: faces 11r, 1rr
RIF(z.subs({ra:1,rb:1,rc:r,d:d11rr^(1/3)}).subs(r=sqrt(2)-1)).lower() # sliding 1 in face 11r
RIF(z.subs({ra:r,rb:1,rc:1,d:d11rr^(1/3)}).subs(r=sqrt(2)-1)).lower() # sliding r in face 11r
RIF(z.subs({ra:1,rb:r,rc:r,d:d11rr^(1/3)}).subs(r=sqrt(2)-1)).lower() # sliding 1 in face 1rr
RIF(z.subs({ra:r,rb:1,rc:r,d:d11rr^(1/3)}).subs(r=sqrt(2)-1)).lower() # sliding r in face 1rr

# tetrahedron of type 1rrr: faces 1rr, rrr
RIF(z.subs({ra:1,rb:r,rc:r,d:d1rrr^(1/3)}).subs(r=sqrt(2)-1)).lower() # sliding 1 in face 1rr
RIF(z.subs({ra:r,rb:1,rc:r,d:d1rrr^(1/3)}).subs(r=sqrt(2)-1)).lower() # sliding r in face 1rr
RIF(z.subs({ra:r,rb:r,rc:r,d:d1rrr^(1/3)}).subs(r=sqrt(2)-1)).lower() # sliding r in face rrr

# tetrahedron of type rrrr: face rrr
RIF(z.subs({ra:r,rb:r,rc:r,d:drrrr^(1/3)}).subs(r=sqrt(2)-1)).lower() # sliding r in face rrr


