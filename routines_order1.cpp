// quite useful for the volume (improves the 3 or 4 decimal)
// the same approach does not seem that useful on Lagrange formula for solid angle.

// return the first partial derivative of the volume taken in (ab,ac,ad,bc,bd,cd)
I diff_volume(I ab, I ac, I ad, I bc, I bd, I cd)
{
    I x0=square(ab);
    I x1=square(ac);
    I x2=square(ad);
    I x3=square(bc);
    I x4=square(bd);
    I x5=square(cd);

    return -ab*(x1*x3+x0*x5-(x1*x4+(x3-x4)*x2-(x0-x1-x2-x3-x4+x5)*x5))/sqrt((x1*x4+(x3-x4)*x2-(x0-x1-x2-x3-x4+x5)*x5)*x0-(x0*x3-(x3-x5)*x2+(x1-x2-x3+x4-x5)*x4)*x1-(x4*x5+(x2+x3-x4-x5)*x2)*x3)/I(12); //50
}

// volume of the block using Mean value theorem around the block center
// CAUTION: return empty if the block center is not tetrahedral, though the block can contains tetrahedra (hence cannot be used to discard non-tetrahedral blocks but only to refine the volume computation) 
I volume_order1(I ab, I ac, I ad, I bc, I bd, I cd)
{
    I ab0=median(ab);
    I ac0=median(ac);
    I ad0=median(ad);
    I bc0=median(bc);
    I bd0=median(bd);
    I cd0=median(cd);

    I V=volume(ab0,ac0,ad0,bc0,bd0,cd0);

    // add partial derivatives, taken on the whole block and multiplied by the distance to the center
    V+=(ab-ab0)*diff_volume(ab,ac,ad,bc,bd,cd);
    V+=(ac-ac0)*diff_volume(ac,ad,ab,cd,bc,bd);
    V+=(ad-ad0)*diff_volume(ad,ab,ac,bd,cd,bc);
    V+=(bc-bc0)*diff_volume(bc,cd,ac,bd,ab,ad);
    V+=(bd-bd0)*diff_volume(bd,ad,cd,ab,bc,ac);
    V+=(cd-cd0)*diff_volume(cd,ac,bc,ad,bd,ab);

    return V;
}

// return the firt partial derivative of lagrange taken in (ab,ac,ad,bc,bd,cd)
I diff_lagrange(I ab, I ac, I ad, I bc, I bd, I cd)
{
    return (square(ab+ac+ad)-square(ab)-square(cd))/I(12);
}

// lagrange of the tetrahedron using Mean value theorem around the block center
I lagrange_order1(I ab, I ac, I ad, I bc, I bd, I cd)
{
    I ab0=median(ab);
    I ac0=median(ac);
    I ad0=median(ad);
    I bc0=median(bc);
    I bd0=median(bd);
    I cd0=median(cd);

    I f=lagrange(ab0,ac0,ad0,bc0,bd0,cd0);

    // sum each partial derivative
    f+=(ab-ab0)*diff_lagrange(ab,ac,ad,bc,bd,cd);
    f+=(ac-ac0)*diff_lagrange(ac,ad,ab,cd,bc,bd);
    f+=(ad-ad0)*diff_lagrange(ad,ab,ac,bd,cd,bc);
    f+=(bc-bc0)*diff_lagrange(bc,cd,ac,bd,ab,ad);
    f+=(bd-bd0)*diff_lagrange(bd,ad,cd,ab,bc,ac);
    f+=(cd-cd0)*diff_lagrange(cd,ac,bc,ad,bd,ab);

    return f;
}


