#include "CCStateSpace.h"
#include "mex.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	if (nrhs != 6)//
	{
		mexErrMsgTxt("the number of inputs is must 6");
	}
	if (nlhs != 2)
	{
		mexErrMsgTxt("the number of outputs is must 2");
	}

	//----------------Inputs-------------------
	double *q1, *q2, *tp,*para, *sp,*step;
	q1= mxGetPr(prhs[0]); //q_s
	q2= mxGetPr(prhs[1]); //q_g
	tp = mxGetPr(prhs[2]); //type
	para = mxGetPr(prhs[3]); //parameters
	sp = mxGetPr(prhs[4]); //space
	step = mxGetPr(prhs[5]); //step

	CCStateSpace space(sp[0], sp[1]);
	CCStateSpace::State q_s(space, q1[0], q1[1], q1[2]);
	CCStateSpace::State q_g(space, q2[0], q2[1], q2[2]);

	
	//---------------End of Inputs--------------
	vector<CCStateSpace::State> vec;
	double err_t=-1.;
	space.interpolate(&err_t, &vec, &q_s, &q_g, *tp, para, *step);

	
	//
	if (!vec.empty())
	{
		size_t n = vec.size();
		plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL); //err
		plhs[1] = mxCreateDoubleMatrix(n, 3, mxREAL); //parameters
		
		
		double *err = mxGetPr(plhs[0]);
		double *state = mxGetPr(plhs[1]);

		*err = err_t;
		for (int i = 0; i < n; ++i)
		{
			state[i] = (vec[i]).x;
			state[n + i] = (vec[i]).y;
			state[2 * n + i] = (vec[i]).theta;
		}
	}
	else
	{
		plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL); //err
		plhs[1] = mxCreateDoubleMatrix(1, 1, mxREAL); //states

		double *err = mxGetPr(plhs[0]);
		double *state = mxGetPr(plhs[1]);
		
		*err = err_t;
		*state = 0;
	}
}