
// return false if x<=y+z (or one of the permutation) does not hold 
bool trig(I x, I y, I z)
{
    return (lower(x)<=upper(y+z) && lower(y)<=upper(x+z) && lower(z)<=upper(x+y));
}

// return false if the facial condition is not satisfied
bool is_facial(I ab, I ac, I ad, I bc, I bd, I cd)
{
    return (trig(ab,ac,bc) && trig(ab,bd,ad) && trig(ac,cd,ad) && trig(bc,bd,cd));
}

// volume of the tetrahedron
I volume(I ab, I ac, I ad, I bc, I bd, I cd)
{
    I x0=square(ab);
    I x1=square(ac);
    I x2=square(ad);
    I x3=square(bc);
    I x4=square(bd);
    I x5=square(cd);

    return sqrt((I(4)*x2*x0-square(x0+x2-x4))*(I(4)*x1*x0-square(x0+x1-x3))-square((x1+x2-x5)*I(2)*x0-(x0+x1-x3)*(x0+x2-x4)))/I(24)/ab; // 21 variables
}

// return f such that tan(A/2)=V/f where V is the volume and A the solid angle in A (Lagrange formula)
I lagrange(I ab, I ac, I ad, I bc, I bd, I cd)
{

    I x0=square(ab);
    I x1=square(ac);
    I x2=square(ad);
    I x3=square(bc);
    I x4=square(bd);
    I x5=square(cd);
    return (I(2)*ab*ac*ad+(x1+x2-x5)*ab+(x0+x2-x4)*ac+(x0+x1-x3)*ad)/I(12);
}

// compute the density using the Lagrange formula for solid angles and V for the volume
I density(I ab, I ac, I ad, I bc, I bd, I cd, I V)
{
    I f,A,d,dt=I(0);

    // we add to dt the contribution of every vertex

    f=lagrange(ab,ac,ad,bc,bd,cd); // denominator in the Lagrange formula
    A=I(2)*atan(V/f); // solid angle in [-pi,pi]
    if (upper(A)<=0) A+=pi_twice<I>(); // solid angle in [0,2pi]
    d=A/V; // usual formula
    if (lower(V/f)>=0) d=intersect(d,hull(I(0),I(2)/f)); // using atan(x)<=x for x>=0 : V simplifies, interesting if V small!
    dt+=d*pow(ra,3)/I(3);

    f=lagrange(bc,bd,ab,cd,ac,ad);
    A=I(2)*atan(V/f);
    if (upper(A)<=0) A+=pi_twice<I>();
    d=A/V;
    if (lower(V/f)>=0) d=intersect(d,hull(I(0),I(2)/f));
    dt+=d*pow(rb,3)/I(3);

    f=lagrange(cd,ac,bc,ad,bd,ab);
    A=I(2)*atan(V/f);
    if (upper(A)<=0) A+=pi_twice<I>();
    d=A/V;
    if (lower(V/f)>=0) d=intersect(d,hull(I(0),I(2)/f));
    dt+=d*pow(rc,3)/I(3);

    f=lagrange(ad,bd,cd,ab,ac,bc);
    A=I(2)*atan(V/f);
    if (upper(A)<=0) A+=pi_twice<I>();
    d=A/V;
    if (lower(V/f)>=0) d=intersect(d,hull(I(0),I(2)/f));
    dt+=d*pow(rd,3)/I(3);

    return dt;
}

// is the block B near an optimal block?
bool is_near_optimal(block B)
{
    // count the number of eps-tight edges
    int nb_contacts=0;
    if (upper(B[0])<lower(ra+rb+eps)) nb_contacts++;
    if (upper(B[1])<lower(ra+rc+eps)) nb_contacts++;
    if (upper(B[2])<lower(ra+rd+eps)) nb_contacts++;
    if (upper(B[3])<lower(rb+rc+eps)) nb_contacts++;
    if (upper(B[4])<lower(rb+rd+eps)) nb_contacts++;
    if (upper(B[5])<lower(rc+rd+eps)) nb_contacts++;

    // the tight uuuu, uurr, urrr
    #if defined(case_1) || defined(case_4)  || defined(case_5) || defined(case_6) || defined(case_7) || defined(case_8)
    if (nb_contacts==6) return true;
    #endif

    // the bi-stretched rrrr
    #if defined(case_9)
    // 4 tight edges and 2 stretched rr-edges
    if (nb_contacts!=4) return false;
    int j=0;
    for(int i=0;i<6;i++)
        if (lower(B[i])>upper(stretched_rrrr_rr-eps)) j++;
    if (j==2) return true;
    #endif

    // the uu-stretched uuur
    #if defined(case_2) || defined(case_3)
    //5 tight edges and 1 stretched uu-edge
    if (nb_contacts!=5) return false;
    for(int i=0;i<6;i++)
        if ((i<2 || i==3) && lower(B[i])>upper(stretched_uuur_uu-eps)) return true; // stretched uu
    #endif

    //in all other cases the block is regular
    return false;
}
