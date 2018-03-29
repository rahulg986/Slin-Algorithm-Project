/*solving 0.5x'Qx + bx subject to: -lambda <= x <= lambda*/
//input Q,b,l,u,start point, iterations
#include "mex.h"
#include "math.h"
#if defined(NAN_EQUALS_ZERO)
#define IsNonZero(d) ((d)!=0.0 || mxIsNaN(d))
#else
#define IsNonZero(d) ((d)!=0.0)
#endif
void mexFunction(int nlhs,mxArray *plhs[],int nrhs, const mxArray *prhs[])
{
    mwSize n,nzmax;
    double *q,*b,*l,*u,*cur_sol,*y;
    double c1,c2,val,change,iteration;
    mwIndex *ir,*jc;
    mwIndex i,iter,nz_col,start,k,row_pos,m;
    nzmax=mxGetNzmax(prhs[0]);
    q=mxGetPr(prhs[0]);
    m=mxGetM(prhs[0]);
    ir=mxGetIr(prhs[0]);
    jc=mxGetJc(prhs[0]);
    b=mxGetPr(prhs[1]);
    l=mxGetPr(prhs[2]);
    u=mxGetPr(prhs[3]);
    cur_sol=mxGetPr(prhs[4]);
    iteration=(int)mxGetScalar(prhs[5]);
    plhs[0]=mxCreateDoubleMatrix(m,1,mxREAL);
    y=mxGetPr(plhs[0]);
    for (i=0;i < m; i++)
    {
       *(y+i)=*(cur_sol+i);
      /* mexprintf("it is:",i,*(y+i));*/
    }
    for (iter=0; iter < iteration; iter++)
    {
        change=0;
        //Traverse columns of sparse matrix Q
        for (i=0; i < m; i++)
        {
            c1=0;
            c2=0;
            nz_col=*(jc+i+1)-*(jc+i);
            start=*(jc+i);
            //Traverse column i
            for (k=0; k < nz_col; k++)
            {
                row_pos=*(ir+start+k);
                if (row_pos==i) c1=*(q+start+k)/2;
                else c2=c2+*(q+start+k)**(y+row_pos);
            }
            c2=c2+*(b+i);
            val=-c2/(2*c1);
            if (val > *(u+i))
            {
                change=change + (*(y+i)-*(u+i))*(*(y+i)-*(u+i));
                *(y+i)=*(u+i);
            }
            else if (val < *(l+i))
            {
                change=change + (*(y+i)-*(l+i))*(*(y+i)-*(l+i));
                *(y+i)=*(l+i);
            }
            else
            {
                change=change + (*(y+i)-val)*(*(y+i)-val);
                *(y+i)=val;
            }
        }
        change=sqrt(change);
        if (change < 1e-3) break;
    }

   
    
    
    
}