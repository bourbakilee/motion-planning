#ifndef CCSTATAESPACE_H
#define CCSTATAESPACE_H
#include "fresnel.h"
#include <boost/math/constants/constants.hpp>
#include <cmath>
#include <vector>
using namespace std;

const double CC_EPS = 1e-8;
//const double CC_CHECK=1e-5;
const double pi = boost::math::constants::pi<double>();
const double half_pi = boost::math::constants::half_pi<double>();
const double two_pi = boost::math::constants::two_pi<double>();
const double root_pi = boost::math::constants::root_pi<double>();
const double root_one_div_pi = boost::math::constants::root_one_div_pi<double>();

class CCStateSpace
{
public:
/*-------------------------State--------------------------------*/
	//checked 10.04.2015
	class State{
	public:
		State() :State(0.,0.,0.){}	//
		State(const CCStateSpace &space) :State(space,0.,0.,0.){}	//
		State(double x_, double y_, double theta_) :x(x_),y(y_),theta(theta_){}	//
		State(const CCStateSpace &space, double x_, double y_, double theta_) :x(x_), y(y_), theta(theta_),  x_o_lp(x_ + space.x_o*cos(theta_) - space.y_o*sin(theta_)), y_o_lp(y_ + space.x_o*sin(theta_) + space.y_o*cos(theta_)), x_o_rp(x_ + space.x_o*cos(theta_) + space.y_o*sin(theta_)), y_o_rp(y_ + space.x_o*sin(theta_) - space.y_o*cos(theta_)), x_o_ln(x_ - space.x_o*cos(theta_) - space.y_o*sin(theta_)), y_o_ln(y_ - space.x_o*sin(theta_) + space.y_o*cos(theta_)), x_o_rn(x_ - space.x_o*cos(theta_) + space.y_o*sin(theta_)), y_o_rn(y_ - space.x_o*sin(theta_) - space.y_o*cos(theta_)){}
		State(const State &q) :x(q.x), y(q.y), theta(q.theta), x_o_lp(q.x_o_lp), y_o_lp(q.y_o_lp), x_o_rp(q.x_o_rp), y_o_rp(q.y_o_rp), x_o_ln(q.x_o_ln), y_o_ln(q.y_o_ln), x_o_rn(q.x_o_rn), y_o_rn(q.y_o_rn){}
		
		void cc_centre(const CCStateSpace &space);
		State local_to_ref(State *global_ref, const CCStateSpace &space);
		State global_by_ref(const State &local_ref);
		//
		double x, y, theta;
		
		//
		double x_o_lp, y_o_lp;
		double x_o_rp, y_o_rp;
		double x_o_ln, y_o_ln;
		double x_o_rn, y_o_rn;
	};
/*------------------------CC-Truns------------------------------*/
	//
	static const int ccTurnType[14][4];
	// verified 11.04.2015
	class CCTurn{
	public:
		CCTurn(const int *type_, const State &start, const CCStateSpace &space):CCTurn(type_, start,2*space.r*sin(space.miu),space){}
		//CCTurn(const int *type_, const State &start, double s, const CCStateSpace &space);	//
		CCTurn(const int *type_, const State &start, double delta_, double sigma_, const CCStateSpace &space);	//topology
		CCTurn(const int *type_, const State &start, double delta_, const CCStateSpace &space); //CCTurn
		CCTurn(const CCTurn &cct);
		CCTurn & operator = (const CCTurn &cct);

		State interpolate(double s, const CCStateSpace &space);	//
		State interpolate(int i, int n, const CCStateSpace &space);	//
		void vec_interpolate(double s, vector<State> *vec, const CCStateSpace &space);
		
		//
		const int *type;	//
		State q_s;
		State q_g;
		double delta;
		double sigma;
		double vlength[3];	//
		double length;
	};
/*------------------------CC-Paths------------------------------*/
	enum CCPathType {
		CC_NON = 0,
		LpSpLp = 1, LnSnLn = 2, RpSpRp = 3, RnSnRn = 4,
		LpSpRp = 5, LnSnRn = 6, RpSpLp = 7, RnSnLn = 8,
		LpRnLp = 9, LnRpLn = 10, RpLnRp = 11, RnLpRn = 12,
		LpRnLn = 13, LnRpLp = 14, RpLnRn = 15, RnLpRp = 16,
		LnRnLp = 17, LpRpLn = 18, RnLnRp = 19, RpLpRn = 20,
		LpRpLnRn = 21, LnRnLpRp = 22, RpLpRnLn = 23, RnLnRpLp = 24,
		LpRnLnRp = 25, LnRpLpRn = 26, RpLnRnLp = 27, RnLpRpLn = 28,
		LpRnSnLn = 29, LnRpSpLp = 30, RpLnSnRn = 31, RnLpSpRp = 32,
		LnSnRnLp = 33, LpSpRpLn = 34, RnSnLnRp = 35, RpSpLpRn = 36,
		LpRnSnRn = 37, LnRpSpRp = 38, RpLnSnLn = 39, RnLpSpLp = 40,
		RnSnRnLp = 41, RpSpRpLn = 42, LnSnLnRp = 43, LpSpLpRn = 44,
		LpRnSnLnRp = 45, LnRpSpLpRn = 46, RpLnSnRnLp = 47, RnLpSpRpLn = 48,
		Topology = 49
	};

	static const int PathType[50][5];
	
	class CCPath{
	public:
		//
		CCPath(State *start, State *goal, CCStateSpace *space);
		//
		void vec_interpolate(double step, vector<State> *vec, const CCStateSpace &space);
		//
		State *q_s;
		State *q_g;	//

		double t, u, v;	//
		vector<CCTurn>cc_path; //
		double cc_length;
		CCPathType cc_type;
		double cc_err;
	};
/*--------------------------------------------------------------*/
	//member functions
	CCStateSpace():CCStateSpace(0.2, 0.04){}
	CCStateSpace(double k, double sigma) :k_max(k), sigma_max(sigma), delta_min(k*k/sigma),theta_i(delta_min/2.), x_i(root_pi / sqrt(sigma_max)*fresnel_c(root_one_div_pi * sqrt(delta_min))), y_i(root_pi / sqrt(sigma_max)*fresnel_s(root_one_div_pi * sqrt(delta_min))), x_o(x_i - sin(theta_i) / k_max), y_o(y_i + cos(theta_i) / k_max), r(sqrt(x_o*x_o + y_o*y_o)), miu(atan(x_o / y_o)){}
	
	//
	CCTurn CC_SP(State *q_s, double t);
	CCTurn CC_SN(State *q_s, double t);
	CCTurn CC_LP(State *q_s, double t);
	CCTurn CC_RP(State *q_s, double t);
	CCTurn CC_LN(State *q_s, double t);
	CCTurn CC_RN(State *q_s, double t);
	//
	void CSC(State *q_s, State *q_g, vector<vector<CCTurn>> *vec_cct, vector<double> *t, vector<double> *u, vector<double> *v, vector<double> *length, vector<double> *err, vector<CCPathType> *type);
	void CCC(State *q_s, State *q_g, vector<vector<CCTurn>> *vec_cct, vector<double> *t, vector<double> *u, vector<double> *v, vector<double> *length, vector<double> *err, vector<CCPathType> *type);
	void CCCC(State *q_s, State *q_g, vector<vector<CCTurn>> *vec_cct, vector<double> *t, vector<double> *u, vector<double> *v, vector<double> *length, vector<double> *err, vector<CCPathType> *type);
	void CCSC(State *q_s, State *q_g, vector<vector<CCTurn>> *vec_cct, vector<double> *t, vector<double> *u, vector<double> *v, vector<double> *length, vector<double> *err, vector<CCPathType> *type);
	void CCSCC(State *q_s, State *q_g, vector<vector<CCTurn>> *vec_cct, vector<double> *t, vector<double> *u, vector<double> *v, vector<double> *length, vector<double> *err, vector<CCPathType> *type);
	void TOPOLOGY(State *q_s, State *q_g, vector<vector<CCTurn>> *vec_cct, vector<double> *t, vector<double> *u, vector<double> *v, vector<double> *length, vector<double> *err, vector<CCPathType> *type);
	//
	bool LPSPLP(State *q_s, State *q_g, vector<CCTurn> *vec_cct, double *t, double *u, double *v, double *err,CCPathType *type);
	bool LNSNLN(State *q_s, State *q_g, vector<CCTurn> *vec_cct, double *t, double *u, double *v, double *err,CCPathType *type);
	bool RPSPRP(State *q_s, State *q_g, vector<CCTurn> *vec_cct, double *t, double *u, double *v, double *err,CCPathType *type);
	bool RNSNRN(State *q_s, State *q_g, vector<CCTurn> *vec_cct, double *t, double *u, double *v, double *err,CCPathType *type);
	bool LPSPRP(State *q_s, State *q_g, vector<CCTurn> *vec_cct, double *t, double *u, double *v, double *err,CCPathType *type);
	bool LNSNRN(State *q_s, State *q_g, vector<CCTurn> *vec_cct, double *t, double *u, double *v, double *err,CCPathType *type);
	bool RPSPLP(State *q_s, State *q_g, vector<CCTurn> *vec_cct, double *t, double *u, double *v, double *err,CCPathType *type);
	bool RNSNLN(State *q_s, State *q_g, vector<CCTurn> *vec_cct, double *t, double *u, double *v, double *err,CCPathType *type);

	bool LPRNLP(State *q_s, State *q_g, vector<CCTurn> *vec_cct, double *t, double *u, double *v, double *err,CCPathType *type);
	bool LNRPLN(State *q_s, State *q_g, vector<CCTurn> *vec_cct, double *t, double *u, double *v, double *err,CCPathType *type);
	bool RPLNRP(State *q_s, State *q_g, vector<CCTurn> *vec_cct, double *t, double *u, double *v, double *err,CCPathType *type);
	bool RNLPRN(State *q_s, State *q_g, vector<CCTurn> *vec_cct, double *t, double *u, double *v, double *err,CCPathType *type);
	bool LPRNLN(State *q_s, State *q_g, vector<CCTurn> *vec_cct, double *t, double *u, double *v, double *err,CCPathType *type);
	bool LNRPLP(State *q_s, State *q_g, vector<CCTurn> *vec_cct, double *t, double *u, double *v, double *err,CCPathType *type);
	bool RPLNRN(State *q_s, State *q_g, vector<CCTurn> *vec_cct, double *t, double *u, double *v, double *err,CCPathType *type);
	bool RNLPRP(State *q_s, State *q_g, vector<CCTurn> *vec_cct, double *t, double *u, double *v, double *err,CCPathType *type);
	bool LNRNLP(State *q_s, State *q_g, vector<CCTurn> *vec_cct, double *t, double *u, double *v, double *err,CCPathType *type);
	bool LPRPLN(State *q_s, State *q_g, vector<CCTurn> *vec_cct, double *t, double *u, double *v, double *err,CCPathType *type);
	bool RNLNRP(State *q_s, State *q_g, vector<CCTurn> *vec_cct, double *t, double *u, double *v, double *err,CCPathType *type);
	bool RPLPRN(State *q_s, State *q_g, vector<CCTurn> *vec_cct, double *t, double *u, double *v, double *err,CCPathType *type);

	bool LPRPLNRN(State *q_s, State *q_g, vector<CCTurn> *vec_cct, double *t, double *u, double *v, double *err,CCPathType *type);
	bool LNRNLPRP(State *q_s, State *q_g, vector<CCTurn> *vec_cct, double *t, double *u, double *v, double *err,CCPathType *type);
	bool RPLPRNLN(State *q_s, State *q_g, vector<CCTurn> *vec_cct, double *t, double *u, double *v, double *err,CCPathType *type);
	bool RNLNRPLP(State *q_s, State *q_g, vector<CCTurn> *vec_cct, double *t, double *u, double *v, double *err,CCPathType *type);
	bool LPRNLNRP(State *q_s, State *q_g, vector<CCTurn> *vec_cct, double *t, double *u, double *v, double *err,CCPathType *type);
	bool LNRPLPRN(State *q_s, State *q_g, vector<CCTurn> *vec_cct, double *t, double *u, double *v, double *err,CCPathType *type);
	bool RPLNRNLP(State *q_s, State *q_g, vector<CCTurn> *vec_cct, double *t, double *u, double *v, double *err,CCPathType *type);
	bool RNLPRPLN(State *q_s, State *q_g, vector<CCTurn> *vec_cct, double *t, double *u, double *v, double *err,CCPathType *type);

	bool LPRNSNLN(State *q_s, State *q_g, vector<CCTurn> *vec_cct, double *t, double *u, double *v, double *err,CCPathType *type);
	bool LNRPSPLP(State *q_s, State *q_g, vector<CCTurn> *vec_cct, double *t, double *u, double *v, double *err,CCPathType *type);
	bool RPLNSNRN(State *q_s, State *q_g, vector<CCTurn> *vec_cct, double *t, double *u, double *v, double *err,CCPathType *type);
	bool RNLPSPRP(State *q_s, State *q_g, vector<CCTurn> *vec_cct, double *t, double *u, double *v, double *err,CCPathType *type);
	bool LNSNRNLP(State *q_s, State *q_g, vector<CCTurn> *vec_cct, double *t, double *u, double *v, double *err,CCPathType *type);
	bool LPSPRPLN(State *q_s, State *q_g, vector<CCTurn> *vec_cct, double *t, double *u, double *v, double *err,CCPathType *type);
	bool RNSNLNRP(State *q_s, State *q_g, vector<CCTurn> *vec_cct, double *t, double *u, double *v, double *err,CCPathType *type);
	bool RPSPLPRN(State *q_s, State *q_g, vector<CCTurn> *vec_cct, double *t, double *u, double *v, double *err,CCPathType *type);
	bool LPRNSNRN(State *q_s, State *q_g, vector<CCTurn> *vec_cct, double *t, double *u, double *v, double *err,CCPathType *type);
	bool LNRPSPRP(State *q_s, State *q_g, vector<CCTurn> *vec_cct, double *t, double *u, double *v, double *err,CCPathType *type);
	bool RPLNSNLN(State *q_s, State *q_g, vector<CCTurn> *vec_cct, double *t, double *u, double *v, double *err,CCPathType *type);
	bool RNLPSPLP(State *q_s, State *q_g, vector<CCTurn> *vec_cct, double *t, double *u, double *v, double *err,CCPathType *type);
	bool RNSNRNLP(State *q_s, State *q_g, vector<CCTurn> *vec_cct, double *t, double *u, double *v, double *err,CCPathType *type);
	bool RPSPRPLN(State *q_s, State *q_g, vector<CCTurn> *vec_cct, double *t, double *u, double *v, double *err,CCPathType *type);
	bool LNSNLNRP(State *q_s, State *q_g, vector<CCTurn> *vec_cct, double *t, double *u, double *v, double *err,CCPathType *type);
	bool LPSPLPRN(State *q_s, State *q_g, vector<CCTurn> *vec_cct, double *t, double *u, double *v, double *err,CCPathType *type);

	bool LPRNSNLNRP(State *q_s, State *q_g, vector<CCTurn> *vec_cct, double *t, double *u, double *v, double *err,CCPathType *type);
	bool LNRPSPLPRN(State *q_s, State *q_g, vector<CCTurn> *vec_cct, double *t, double *u, double *v, double *err,CCPathType *type);
	bool RPLNSNRNLP(State *q_s, State *q_g, vector<CCTurn> *vec_cct, double *t, double *u, double *v, double *err,CCPathType *type);
	bool RNLPSPRPLN(State *q_s, State *q_g, vector<CCTurn> *vec_cct, double *t, double *u, double *v, double *err,CCPathType *type);

	
	//
	static double distance(const State &q1, const State &q2,const CCStateSpace& space);
	static double distance(State *q1, State *q2,const CCStateSpace& space);
	double distance(State *q1, State *q2);
	bool interpolate(double *err, vector<State> *vec, State *q_s, State *q_g, int type, double *parameters, double step);
	// member datas
	const double k_max;
	const double sigma_max;
	//
	const double delta_min;
	const double theta_i,x_i, y_i;
	const double x_o, y_o, r, miu;
	
};

#endif