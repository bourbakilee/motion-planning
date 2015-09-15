#include "CCStateSpace.h"
#include "mex.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	if (nrhs != 4)//
	{
		mexErrMsgTxt("the number of inputs is must 4");
	}
	if (nlhs != 5)
	{
		mexErrMsgTxt("the number of outputs is must 5");
	}

	//----------------Inputs-------------------
	double *q1, *q2, *s, *sp;
	q1= mxGetPr(prhs[0]); //q_s
	q2= mxGetPr(prhs[1]); //q_g
	s = mxGetPr(prhs[2]); //step
	sp = mxGetPr(prhs[3]); //space

	CCStateSpace space(sp[0],sp[1]);
	CCStateSpace::State q_s(space, q1[0], q1[1], q1[2]);
	CCStateSpace::State q_g(space, q2[0], q2[1], q2[2]);
	double step = s[0];
	//---------------End of Inputs--------------

	CCStateSpace::CCPath path(&q_s, &q_g,&space);
	vector<CCStateSpace::State> vec;
	path.vec_interpolate(step, &vec, space);

	
	//
	if (!vec.empty())
	{
		size_t n = vec.size();
		plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL); //type
		plhs[1] = mxCreateDoubleMatrix(1, 3, mxREAL); //parameters
		plhs[2] = mxCreateDoubleMatrix(1, 1, mxREAL); //length
		plhs[3] = mxCreateDoubleMatrix(1, 1, mxREAL); //err
		plhs[4] = mxCreateDoubleMatrix(n, 3, mxREAL); //states
		double *type = mxGetPr(plhs[0]);
		double *parameters = mxGetPr(plhs[1]);
		double *length = mxGetPr(plhs[2]);
		double *err = mxGetPr(plhs[3]);
		double *state = mxGetPr(plhs[4]);

		*type = path.cc_type;
		parameters[0] = path.t;
		parameters[1] = path.u;
		parameters[2] = path.v;
		*length = path.cc_length;
		*err = path.cc_err;
		for (int i = 0; i < n; ++i)
		{
			state[i] = (vec[i]).x;
			state[n + i] = (vec[i]).y;
			state[2 * n + i] = (vec[i]).theta;
		}
	}
	else
	{
		plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL); //type
		plhs[1] = mxCreateDoubleMatrix(1, 3, mxREAL); //parameters
		plhs[2] = mxCreateDoubleMatrix(1, 1, mxREAL); //length
		plhs[3] = mxCreateDoubleMatrix(1, 1, mxREAL); //err
		plhs[4] = mxCreateDoubleMatrix(1, 1, mxREAL); //states
		double *type = mxGetPr(plhs[0]);
		double *parameters = mxGetPr(plhs[1]);
		double *length = mxGetPr(plhs[2]);
		double *err = mxGetPr(plhs[3]);
		double *state = mxGetPr(plhs[4]);
		*type = path.cc_type;
		parameters[0] = path.t;
		parameters[1] = path.u;
		parameters[2] = path.v;
		*length = path.cc_length;
		*err = path.cc_err;
		*state = 0;
	}
}