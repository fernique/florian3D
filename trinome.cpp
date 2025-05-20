// INPUT: three intervals A,B,C ; optional intervals D and X
// OUTPUT: a vector of 0 to 2 intervals containing the roots of AX^2+BX+C
// OPTIONAL INPUT:
// interval D used for the discriminant of AX^2+BX+C (instead of computing it).
// interval X: output only the roots that intersect X (X=[0,r] for R, X=[ra+rc,ra+rc+2*r]^2) for ac^2)

auto trinome(I A, I B, I C, I D=I(1,0),I X=I(0,INFINITY))
{
    std::vector<I> Xi; // list of valid roots (0 to 2)

    // no chance for AX^2+BX+C to have a root in the wanted interval -> empty interval
    if (!zero_in(A*square(X)+B*X+C)) return Xi;

    // the discriminant must be nonnegative
    if (empty(D)) D=square(B)-I(4)*A*C; // compute the discriminant if not given as an optional parameter
    if (upper(D)<0) return Xi;
    D=intersect(D,I(0,INFINITY));

    // usual formula (useless if 0 in A)
    I X1=(-B+sqrt(D))/I(2)/A;
    I X2=(-B-sqrt(D))/I(2)/A;

    // alternative formula if 0 in A (can also be useful for small A)
    I x=I(4)*A*C/square(B); // no tetrahedron with A=B=0 => 0 cannot be in both A and B forever
    x=intersect(x,I(-INFINITY,1)); // x<=1 iff D>=0
    // the first root tends towards -C/B
    X1=intersect(X1,-C/B*I(2)/(I(1)+sqrt(I(1)-x)));
    // the second one to infinity
    I X2_inv=-A/B*I(2)/(I(1)+sqrt(I(1)-x)); // interval [a,b] which contains 1/X2
    float b_inv=lower(I(1)/I(upper(X2_inv))); // lower bound on 1/b
    X2=intersect(X2,I(b_inv,INFINITY)); // X2 (if positive) is at least 1/b

    if (overlap(X1,X2))  // both roots cannot be distinguished
    {
        X1=hull(X1,X2);
        if (overlap(X1,X)) Xi.push_back(X1);
//        Xi.push_back(hull(X1,X2));
    }
    else
    {
        if (overlap(X1,X)) Xi.push_back(X1);
        if (overlap(X2,X)) Xi.push_back(X2);
    }
    return Xi;
}
