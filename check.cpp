// select case_1 to case_9 before compiling (g++ check.cpp -o check -O3)
// adjust the constant maxThreads to the number of cores of the used computer
#define case_1

/*************************************************************/
// Interval arithmetic

#include <boost/numeric/interval.hpp>
using namespace boost::numeric;
using namespace interval_lib;
// Fix rounding policies for the transcendental functions
typedef interval<double, policies<save_state<rounded_transc_std<double>>, checking_base<double>>> I;

/*************************************************************/
// Multithreading

#include <boost/asio/post.hpp>
#include <boost/asio/thread_pool.hpp>

const int maxThreads = 8; // max nb of threads in the pool (8 for my 8-core Dell laptop)
boost::asio::thread_pool pool(maxThreads); // pool of threads

const int maxWaitingThreads=maxThreads; // max nb threads in the queue to enter the pool
std::atomic_int nWaitingThreads=0; 

/*************************************************************/

#include <array>
#include <vector>
#include <math.h>

// block: edge lengths ab,ac,ad,bc,bd,cd and radius R of the support sphere
typedef std::array<I,7> block;

// radii of spheres
const I r=sqrt(I(2))-I(1);
const I u=I(1);

// maximal recursion depth
const int maxDepth=100;
std::atomic_int maxReachedDepth=0; 

// for statistics
unsigned long nBlocks=0;
std::array<std::array<unsigned long,3>,maxDepth> stats;
std::mutex stats_mutex;

/*************************************************************/
// case dependent constants (radii of the sphere, density, eps...)

#if defined(case_1) // uuuu
const I ra=u, rb=u, rc=u, rd=u; // radii of the spheres
const I delta_bound=I(779635570044252)/I(1000000000000000); // lower bound on the claimed maximal density
const I eps=I(1)/I(80); // neighborhood of the optimal

#elif defined(case_2) || defined(case_3) // uuur
const I ra=u, rb=u, rc=u, rd=r;
const I delta_bound=I(812542027810834)/I(1000000000000000);
const I eps=I(1)/I(499);
const I stretched_uuur_uu=I(4)*sqrt(I(2)*r/(I(1)+I(2)*r)); // length 2.69245416472622 of the stretched uu edge

#elif defined(case_4) || defined(case_5) || defined(case_6) // uurr
const I ra=u, rb=u, rc=r, rd=r;
const I delta_bound=I(810466032832072)/I(1000000000000000);
const I eps=I(1)/I(317);

#elif defined(case_7) || defined(case_8) // urrr
const I ra=u, rb=r, rc=r, rd=r;
const I delta_bound=I(806503318194810)/I(1000000000000000);
const I eps=I(1)/I(410);

#elif defined(case_9) // rrrr
const I ra=r, rb=r, rc=r, rd=r;
const I delta_bound=I(784688454045207)/I(1000000000000000);
const I eps=I(1)/I(243);
const I stretched_rrrr_rr=r*sqrt(I(2)*sqrt(I(6))+I(6)); // length 1.3674681889063 of the two stretched rr edges
#endif

/*************************************************************/
// Auxiliary files

// routines for tetrahedra (volume, angles, density etc)
#include "routines.cpp"
// more involved routine (in particular for the volume)
#include "routines_order1.cpp"
// to compute the root of a quadratic polynomial
#include "trinome.cpp"
// to update R in a newly created block
#include "radius.cpp"

/*************************************************************/
// return the type of the block B:
// 0=near optimal, 1=invalid, 2=hollow
// -1 : the block B has to be further divided in order to determine its type

int type(block B)
{    
    // block near an optimal block?
    if (is_near_optimal(B)) return 0;

    // the support sphere is too large for a FM-tetrahedron
    // actually useless because radius returns only radii that intersect [0,r])
//    if (lower(B[6])>=upper(r)) return 1;

    // compute the volume V of the block (usual formula)
    I V=volume(B[0],B[1],B[2],B[3],B[4],B[5]);

    // Refine the volume using order 1 formula
    // It is usually better for small blocks
    // however it yields an empty interval if the center of the block is not tetrahedral
    // even if the block itself can contain tetrahedra -> just ignore it in such a case
    I V1=volume_order1(B[0],B[1],B[2],B[3],B[4],B[5]);
    if (!empty(V1)) V=intersect(V,V1);

    // non-tetrahedral block?
    if (!is_facial(B[0],B[1],B[2],B[3],B[4],B[5]) || empty(V)) return 1;
    
    // compute the density d of the block
    I d=density(B[0],B[1],B[2],B[3],B[4],B[5],V);
    
    // hollow block? (density lower than the claimed bound)
    if (upper(d)<=lower(delta_bound)) return 2;

    // if it is not clear whether R is in [0,r] or not (i.e. the interval for R contains 0 or r)
    // then the block must be further subdivided in smaller blocks whose type can be determined
    if (zero_in(B[6]) || zero_in(B[6]-r)) return -1;
    
    // at this point, we are sure that B is a block of FM-tetrahedra.
    // if the density is higher than the claimed bound, then we get a counter-example to the main theorem.
    // Spoiler: this never happens!
    if (lower(d)>upper(delta_bound))
    {
        printf("Found a block with density at least %f:\n",lower(d));
        throw "";
    }
    // otherwise it means we cannot yet compare the density of B with the claimed maximal density
    // the block has to be further subdivided
    return -1;
}


/*************************************************************/
// halve the block B along its widest side
auto halve_block(block B)
{
    int i=0;
    for(int j=1;j<6;j++)
        if (width(B[j])>width(B[i])) i=j;
    std::pair<I,I> e=bisect(B[i]);
    std::vector<block> Bi;
    B[i]=e.first;
    Bi.push_back(B);
    B[i]=e.second;
    Bi.push_back(B);
    return Bi;
}

/*************************************************************/
// recursively subdivide a block whose type cannot yet be determined

void explore(block B, int depth=0);

// auxiliary function in order to count (and bound) the nb of waiting threads
void aux_explore(block B, int depth)
{
    nWaitingThreads--;
    explore(B, depth);
}

void explore(block B, int depth)
{
    int t=type(B);

    if (depth>maxReachedDepth)
        maxReachedDepth=depth;
 
    if (t>=0) // if the block can be ruled out it's great!
    {
        stats_mutex.lock();
        nBlocks++;
        stats[depth][t]+=1;
        if (nBlocks%1000000==0) printf("%lu millions of blocks checked\n",nBlocks/1000000);
        stats_mutex.unlock();
    }
    else // otherwise split it into two smaller blocks and continue on each of them
        for (auto& B_halved: halve_block(B))
            for (auto& B_halved_updated: update_block(B_halved))
                if (nWaitingThreads>maxWaitingThreads) // too much threads are waiting -> direct recursive call
                {
                    explore(B_halved_updated, depth+1);
                }
                else // not that much waiting threads -> post the recursive call in the thread pool
                {
                    nWaitingThreads++;
                    boost::asio::post(pool, [B_halved_updated, depth]{aux_explore(B_halved_updated, depth+1); });
                }
}

/*************************************************************/
// initialize the first block (case dependent)

block initial_block()
{
    block B;
    // first default edge length
    I gap=hull(I(0),I(2)*r);
    B[0]=ra+rb+gap;
    B[1]=ra+rc+gap;
    B[2]=ra+rd+gap;
    B[3]=rb+rc+gap;
    B[4]=rb+rd+gap;
    B[5]=rc+rd+gap;
    // then update edge length depending on which contacts are assumed
    #if defined(case_1) || defined(case_2) || defined(case_4)
    B[0]=ra+rb; // contact along ab
    #elif defined(case_3) || defined(case_5) || defined(case_7)
    B[2]=ra+rd; // contact along ad
    #elif defined(case_6) || defined(case_8) || defined(case_9)
    B[5]=rc+rd; // contact along cd
    #endif
    return B;
}

/*************************************************************/
int main(int argc, char *argv[])
{
    // define the initial block B
    block B=initial_block();

    // recursively subdivide B until the type of every sub-blocks has been determined
    explore(B);
    pool.join(); // wait for all the threads to be terminated

    // print stats
    printf("Number of cuts and corresponding near-optimal, invalid & hollow blocks discarded:\n");
    for(int i;i<maxReachedDepth;i++)
    {
        printf("%d, ",i);
        for(int j=0;j<3;j++) printf("%lu, ",stats[i][j]);
        printf("\n");
    }

    printf("%lu blocks checked\n",nBlocks);

    return 0;
}
