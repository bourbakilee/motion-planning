/*
LI Yunsheng, 2015
*/

#include "CCStateSpace.h"

inline double mod2pi(double x) //-pi~pi
{
    double v = fmod(x,two_pi);
    if (v < -pi)
        v += two_pi;
    else if (v > pi)
        v -= two_pi;
    return v;
}

inline double mod2pi2(double x) //0~2pi
{
    double v = fmod(x,two_pi);
    if (v < 0)
        v += two_pi;
    return v;
}

inline double mod2pi3(double x) //-2pi~0
{
    double v = fmod(x, two_pi);
    if (v > 0)
        v -= two_pi;
    return v;
}

inline double D(double u)
{
    return fresnel_c(sqrt(2 * u)*root_one_div_pi)*cos(u) + fresnel_s(sqrt(2 * u)*root_one_div_pi)*sin(u);
}

inline double d1(double x, const CCStateSpace &space) //x>= space.theta_i
{
    return 2/space.k_max*sqrt(2*pi*x)*D(x)*sin(x)/cos(2*x);
}

inline double d2(double x, const CCStateSpace &space) //0<=x<=space.theta_i
{
    return 2*sqrt(pi/space.sigma_max)*D(x)*sin(x)/cos(2*x);
}

double solve_alpha(double d, double d_m, const CCStateSpace &space) //
{
    double alpha=space.theta_i;
    double d_t;
    
    if(d>d_m+CC_EPS)
    {
        double a=space.theta_i, b=pi/4-CC_EPS;
        do
        {
            alpha=(a+b)/2;
            d_t=d1(alpha, space);
            if(d_t>d)
                b=alpha;
            else
                a=alpha;
        } while (abs(d_t - d) > CC_EPS);
    }
    else if(d<d_m-CC_EPS)
    {
        double a=0, b= space.theta_i;
        do
        {
            alpha=(a+b)/2;
            d_t=d2(alpha, space);
            if(d_t>d)
                b=alpha;
            else
                a=alpha;
        } while (abs(d_t - d) > CC_EPS);
    }
    return alpha;
}

inline double dist2d(double dy,double dx)
{
    return sqrt(dx*dx+dy*dy);
}

inline double dist2(double dy,double dx)
{
    return dy*dy+dx*dx;
}

double CCStateSpace::distance(const State &q1, const State &q2, const CCStateSpace &space)
{
    return dist2d(q2.y-q1.y,q2.x-q1.x)+abs(mod2pi(q2.theta-q1.theta))/space.k_max;
}

double CCStateSpace::distance(State *q1, State *q2,const CCStateSpace& space)
{
    return dist2d(q2->y-q1->y,q2->x-q1->x)+abs(mod2pi(q2->theta-q1->theta))/space.k_max;
}

double CCStateSpace::distance(State *q1, State *q2)
{
    return CCStateSpace::distance(q1,q2,*this);
}

//checked 10.04.2015, {Forward/Backward, Left/Right, Head: Increase/Decrease, Default/Reverse at middle circle/No middle circle}
const int CCStateSpace::ccTurnType[14][4] = {
    {1,0,0,0},  //SP                0
    {-1,0,0,0}, //SN                1
    {1,1,1,1},  //LP                2
    {1,1,1,-1}, //                  3
    {1,1,1,0},  //Lp                4
    {1,-1,-1,1},    //RP            5
    {1,-1,-1,-1},   //              6
    {1,-1,-1,0},    //              7
    {-1,1,-1,1},    //LN            8
    {-1,1,-1,-1},   //              9
    {-1,1,-1,0},    //              10
    {-1,-1,1,1},    //RN            11
    {-1,-1,1,-1},   //              12
    {-1,-1,1,0} //                  13
};

//Lp-1, Rp-2, Ln-3, Rn-4, Sp-5, Sn-6, NON-0
const int CCStateSpace::PathType[50][5]={
    {0,0,0,0,0},    //0
    {1,5,1,0,0},    //1
    {3,6,3,0,0},    //2
    {2,5,2,0,0},    //3
    {4,6,4,0,0},    //4
    {1,5,2,0,0},    //5
    {3,6,4,0,0},    //6
    {2,5,1,0,0},    //7
    {4,6,3,0,0},    //8
    {1,4,1,0,0},    //9
    {3,2,3,0,0},    //10
    {2,3,2,0,0},    //11
    {4,1,4,0,0},    //12
    {1,4,3,0,0},    //13
    {3,2,1,0,0},    //14
    {2,3,4,0,0},    //15
    {4,1,2,0,0},    //16
    {3,4,1,0,0},    //17
    {1,2,3,0,0},    //18
    {4,3,2,0,0},    //19
    {2,1,4,0,0},    //20
    {1,2,3,4,0},    //21
    {3,4,1,2,0},    //22
    {2,1,4,3,0},    //23
    {4,3,2,1,0},    //24
    {1,4,3,2,0},    //25
    {3,2,1,4,0},    //26
    {2,3,4,1,0},    //27
    {4,1,2,3,0},    //28
    {1,4,6,3,0},    //29
    {3,2,5,1,0},    //30
    {2,3,6,4,0},    //31
    {4,1,5,2,0},    //32
    {3,6,4,1,0},    //33
    {1,5,2,3,0},    //34
    {4,6,3,2,0},    //35
    {2,5,1,4,0},    //36
    {1,4,6,4,0},    //37
    {3,2,5,2,0},    //38
    {2,3,6,3,0},    //39
    {4,1,5,1,0},    //40
    {4,6,4,1,0},    //41
    {2,5,2,3,0},    //42
    {3,6,3,2,0},    //43
    {1,5,1,4,0},    //44
    {1,4,6,3,2},    //45
    {3,2,5,1,4},    //46
    {2,3,6,4,1},    //47
    {4,1,5,2,3},    //48
    {7,7,7,7,7}     //49
};

//checked 10.04.2015
void CCStateSpace::State::cc_centre(const CCStateSpace &space)
{
    x_o_lp = x + space.x_o*cos(theta) - space.y_o*sin(theta);
    y_o_lp = y + space.x_o*sin(theta) + space.y_o*cos(theta);
    x_o_rp = x + space.x_o*cos(theta) + space.y_o*sin(theta);
    y_o_rp = y + space.x_o*sin(theta) - space.y_o*cos(theta);
    x_o_ln = x - space.x_o*cos(theta) - space.y_o*sin(theta);
    y_o_ln = y - space.x_o*sin(theta) + space.y_o*cos(theta);
    x_o_rn = x - space.x_o*cos(theta) + space.y_o*sin(theta);
    y_o_rn = y - space.x_o*sin(theta) - space.y_o*cos(theta);
}

//the relative coordinate of q_s to q_g, used to calculate topological path, checked 10.04.2015
CCStateSpace::State  CCStateSpace::State::local_to_ref(State *global_ref, const CCStateSpace &space)
{
    double theta_r = mod2pi(global_ref->theta - theta); //-pi~pi
    double c = cos(theta), s = sin(theta);
    double x_r = (global_ref->x - x)*c + (global_ref->y - y)*s;
    double y_r = -(global_ref->x - x)*s + (global_ref->y - y)*c;
    return CCStateSpace::State(space, x_r, y_r, theta_r);   //
}

//the global coordinate, used to interpolation, checked 10.04.2015
CCStateSpace::State  CCStateSpace::State::global_by_ref(const State &local_ref)
{
    double theta_g = mod2pi2(theta+local_ref.theta);
    double c = cos(theta), s = sin(theta);
    double x_g = x + local_ref.x*c - local_ref.y*s;
    double y_g = y + local_ref.x*s + local_ref.y*c;
    return CCStateSpace::State(x_g, y_g, theta_g);  //
}


//checked 10.04.2015, verified 11.04.2015
//0<delta<delta_min; or topological path(0<delta<pi/2)
//topological path, sigma={double d = D(delta / 2); double s = sin(delta/2+space.miu); sigma = pi*d*d / (space.r*space.r*s*s);}
CCStateSpace::CCTurn::CCTurn(const int *type_, const State &start, double delta_, double sigma_, const CCStateSpace &space) :type(type_),q_s(start),delta(delta_),sigma(sigma_)
{
    double c1 = cos(delta / 2), s1 = sin(delta / 2);
    double d = 2 * root_pi / sqrt(sigma)*D(delta/2);
    double theta_g = mod2pi2(q_s.theta + type[2]*delta);
    double c2 = cos(q_s.theta), s2 = sin(q_s.theta);
    
    vlength[0] = vlength[2] = sqrt(delta / sigma);
    vlength[1] = 0.;
    length = vlength[0] + vlength[2];
    q_g = CCStateSpace::State(space, q_s.x + type[0] * c2*d*c1 - type[1] * s2*d*s1, q_s.y + type[0] * s2*d*c1 + type[1] * c2*d*s1, theta_g);
}

//checked 10.04.2015, verified 11.04.2015
//delta==0; or delta_min<delta<delta_min+pi; or delta_min+pi <delta <2*pi
CCStateSpace::CCTurn::CCTurn(const int * type_, const State &start, double delta_, const CCStateSpace &space) :type(type_), q_s(start), delta(delta_), sigma(space.sigma_max)
{
    //Line,if delta<CC_EPS, CCTurn(q_s,2*r*sin(miu),space), verified 11.04.2015
    if (type[1] == 0)
    {
        q_g = CCStateSpace::State(space, start.x + type[0] * delta*cos(start.theta), start.y + type[0] * delta*sin(start.theta), start.theta);
        vlength[0] = delta; //delta >=0
        vlength[1] = vlength[2] = 0.;
        length = vlength[0];
    }
    else //CCTurn
    {
        double theta_g = mod2pi2(q_s.theta + type[2] * delta);
        double theta_t;
        vlength[0] = vlength[2] = space.k_max / space.sigma_max;
        vlength[1] = space.delta_min <= delta && delta < space.delta_min + pi ? (delta - space.delta_min) / space.k_max : (two_pi - delta + space.delta_min) / space.k_max;
        length = vlength[0] + vlength[1] + vlength[2];

        //
        if (type[0] == 1 && type[1] == 1)   //LP
        {
            theta_t = mod2pi2(theta_g + space.miu - half_pi);   //
            q_g = CCStateSpace::State(space, q_s.x_o_lp + space.r*cos(theta_t), q_s.y_o_lp + space.r*sin(theta_t), theta_g);
        }
        else if (type[0] == 1 && type[1] == -1) //RP
        {
            theta_t = mod2pi2(theta_g - space.miu + half_pi);   //
            q_g = CCStateSpace::State(space, q_s.x_o_rp + space.r*cos(theta_t), q_s.y_o_rp + space.r*sin(theta_t), theta_g);
        }
        else if (type[0] == -1 && type[1] == 1) //LN
        {
            theta_t = mod2pi2(theta_g - space.miu - half_pi);   //
            q_g = CCStateSpace::State(space, q_s.x_o_ln + space.r*cos(theta_t), q_s.y_o_ln + space.r*sin(theta_t), theta_g);
        }
        else    //RN
        {
            theta_t = mod2pi2(theta_g + space.miu + half_pi);   //
            q_g = CCStateSpace::State(space, q_s.x_o_rn + space.r*cos(theta_t), q_s.y_o_rn + space.r*sin(theta_t), theta_g);
        }
    }
}

//verified 11.04.2015
CCStateSpace::CCTurn::CCTurn(const CCTurn &cct) :type(cct.type), q_s(cct.q_s), q_g(cct.q_g), delta(cct.delta), sigma(cct.sigma), length(cct.length)
{
    vlength[0] = cct.vlength[0];
    vlength[1] = cct.vlength[1];
    vlength[2] = cct.vlength[2];
}
CCStateSpace::CCTurn & CCStateSpace::CCTurn::operator = (const CCTurn &cct)
{
    type = cct.type;
    q_s = cct.q_s;
    q_g = cct.q_g;
    delta = cct.delta;
    sigma = cct.sigma;
    length = cct.length;
    vlength[0] = cct.vlength[0];
    vlength[1] = cct.vlength[1];
    vlength[2] = cct.vlength[2];
    
    return *this;
}

//checked 10.04.2015, verified 11.04.2015
CCStateSpace::State  CCStateSpace::CCTurn::interpolate(double s, const CCStateSpace &space) //0<=s<=length
{
    if (s < CC_EPS) return q_s;
    else if (abs(length - s) < CC_EPS) return q_g;
    else if (CC_EPS <= s && s <= length - CC_EPS)
    {
        if (type[1] == 0)   //Line, checked 10.04.2015
            return CCStateSpace::State(q_s.x + type[0] * s*cos(q_s.theta), q_s.y + type[0] * s*sin(q_s.theta), q_s.theta);
        else{
            double tmp = root_pi / sqrt(sigma); //
            if (s <= vlength[0])
                return q_s.global_by_ref(CCStateSpace::State(type[0] * tmp*fresnel_c(s / tmp), type[1] * tmp*fresnel_s(s / tmp), type[2] * sigma*s*s / 2));
            else if (type[3] != 0 && s <= vlength[0] + vlength[1]) //
            {
                double theta_n = mod2pi2(q_s.theta + type[2] * space.theta_i + type[2]*type[3] * (s - vlength[0])*space.k_max); //modified 10.04.2015
                if (type[0] == 1 && type[1] == 1)   //LP
                    return CCStateSpace::State(q_s.x_o_lp + cos(theta_n - half_pi) / space.k_max, q_s.y_o_lp + sin(theta_n - half_pi) / space.k_max, theta_n);
                else if (type[0] == 1 && type[1] == -1) //RP
                    return CCStateSpace::State(q_s.x_o_rp + cos(theta_n + half_pi) / space.k_max, q_s.y_o_rp + sin(theta_n + half_pi) / space.k_max, theta_n);
                else if (type[0] == -1 && type[1] == 1) //LN
                    return CCStateSpace::State(q_s.x_o_ln + cos(theta_n - half_pi) / space.k_max, q_s.y_o_ln + sin(theta_n - half_pi) / space.k_max, theta_n);
                else    //RN
                    return CCStateSpace::State(q_s.x_o_rn + cos(theta_n + half_pi) / space.k_max, q_s.y_o_rn + sin(theta_n + half_pi) / space.k_max, theta_n);
            }
            else
                return q_g.global_by_ref(CCStateSpace::State(-type[0] * tmp*fresnel_c((length - s) / tmp), type[1] * tmp*fresnel_s((length - s) / tmp), -type[2] * sigma*(length - s)*(length - s) / 2));
        }
    }
    else
        return CCStateSpace::State(); //added 11.04.2015
}

//verified 11.04.2015
CCStateSpace::State  CCStateSpace::CCTurn::interpolate(int i, int n, const CCStateSpace &space)
{
    if (i == 0)
        return q_s;
    else if (i == n)
        return q_g;
    else
        return this->interpolate(length*i / n, space);
}

//verified 11.04.2015
void CCStateSpace::CCTurn::vec_interpolate(double s, vector<CCStateSpace::State> *vec, const CCStateSpace &space) //s>>CC_EPS
{
    //s>CC_EPS --step
    if(s>CC_EPS)
    {
        double s_s=0.;
        if(vec->empty()) vec->push_back(q_s);
        if(vlength[0]>CC_EPS)
        {
            for(int i=1;i<ceil(vlength[0]/s);++i)
            {
                s_s=i*s;
                vec->push_back(this->interpolate(s_s,space));
            }
            if(s_s<vlength[0]-CC_EPS) vec->push_back(this->interpolate(vlength[0],space));
        }
        if(vlength[1]>CC_EPS)
        {
            for(int i=1;i<ceil(vlength[1]/s);++i)
            {
                s_s=vlength[0]+i*s;
                vec->push_back(this->interpolate(s_s,space));
            }
            if(s_s<vlength[0]+vlength[1]-CC_EPS) vec->push_back(this->interpolate(vlength[0]+vlength[1],space));
        }
        if(vlength[2]>CC_EPS)
        {
            for(int i=1;i<ceil(vlength[2]/s);++i)
            {
                s_s=vlength[0]+vlength[1]+i*s;
                vec->push_back(this->interpolate(s_s,space));
            }
        if(s_s<length-CC_EPS) vec->push_back(q_g);
        }
    }
    
}


//-----------------------------------Primitive CCTurn-----------------------------------------
//added 11.04.2015
//0<=t<=2*pi
CCStateSpace::CCTurn CCStateSpace::CC_LP(CCStateSpace::State *q_s, double t)
{
    if(t < CC_EPS || t > two_pi - CC_EPS)
    {
        return CCStateSpace::CCTurn(CCStateSpace::ccTurnType[0], *q_s, *this); //line
    }
    else if (t >= CC_EPS && t <= this->delta_min)
    {
        double tmp = D(t / 2) / (this->r*sin(t / 2 + this->miu));
        double sigma_t=pi*tmp*tmp;
        return CCStateSpace::CCTurn(CCStateSpace::ccTurnType[4], *q_s, t, sigma_t, *this); //two clothoids
    }
    else if (t > this->delta_min && t <= this->delta_min + pi)
    {
        return CCStateSpace::CCTurn(CCStateSpace::ccTurnType[2], *q_s, t, *this); //two clothoids and arc
    }
    else
    {
        return CCStateSpace::CCTurn(CCStateSpace::ccTurnType[3], *q_s, t, *this); //two clothoids and arc
    }
}
//added 11.04.2015
CCStateSpace::CCTurn CCStateSpace::CC_RP(CCStateSpace::State *q_s, double t)
{
    if(t < CC_EPS || t > two_pi - CC_EPS)
    {
        return CCStateSpace::CCTurn(CCStateSpace::ccTurnType[0], *q_s, *this); //line
    }
    else if (t >= CC_EPS && t <= this->delta_min)
    {
        double tmp = D(t / 2) / (this->r*sin(t / 2 + this->miu));
        double sigma_t=pi*tmp*tmp;
        return CCStateSpace::CCTurn(CCStateSpace::ccTurnType[7], *q_s, t, sigma_t, *this); //two clothoids
    }
    else if (t > this->delta_min && t <= this->delta_min + pi)
    {
        return CCStateSpace::CCTurn(CCStateSpace::ccTurnType[5], *q_s, t, *this); //two clothoids and arc
    }
    else
    {
        return CCStateSpace::CCTurn(CCStateSpace::ccTurnType[6], *q_s, t, *this); //two clothoids and arc
    }
}
//added 11.04.2015
CCStateSpace::CCTurn CCStateSpace::CC_LN(CCStateSpace::State *q_s, double t)
{
    if(t < CC_EPS || t > two_pi - CC_EPS)
    {
        return CCStateSpace::CCTurn(CCStateSpace::ccTurnType[1], *q_s, *this); //line
    }
    else if (t >= CC_EPS && t <= this->delta_min)
    {
        double tmp = D(t / 2) / (this->r*sin(t / 2 + this->miu));
        double sigma_t=pi*tmp*tmp;
        return CCStateSpace::CCTurn(CCStateSpace::ccTurnType[10], *q_s, t, sigma_t, *this); //two clothoids
    }
    else if (t > this->delta_min && t <= this->delta_min + pi)
    {
        return CCStateSpace::CCTurn(CCStateSpace::ccTurnType[8], *q_s, t, *this); //two clothoids and arc
    }
    else
    {
        return CCStateSpace::CCTurn(CCStateSpace::ccTurnType[9], *q_s, t, *this); //two clothoids and arc
    }
}
//added 11.04.2015
CCStateSpace::CCTurn CCStateSpace::CC_RN(CCStateSpace::State *q_s, double t)
{
    if(t < CC_EPS || t > two_pi - CC_EPS)
    {
        return CCStateSpace::CCTurn(CCStateSpace::ccTurnType[1], *q_s, *this); //line
    }
    else if (t >= CC_EPS && t <= this->delta_min)
    {
        double tmp = D(t / 2) / (this->r*sin(t / 2 + this->miu));
        double sigma_t=pi*tmp*tmp;
        return CCStateSpace::CCTurn(CCStateSpace::ccTurnType[13], *q_s, t, sigma_t, *this); //two clothoids
    }
    else if (t > this->delta_min && t <= this->delta_min + pi)
    {
        return CCStateSpace::CCTurn(CCStateSpace::ccTurnType[11], *q_s, t, *this); //two clothoids and arc
    }
    else
    {
        return CCStateSpace::CCTurn(CCStateSpace::ccTurnType[12], *q_s, t, *this); //two clothoids and arc
    }
}
//added 11.04.2015
CCStateSpace::CCTurn CCStateSpace::CC_SP(CCStateSpace::State *q_s, double t)
{
    return CCStateSpace::CCTurn(ccTurnType[0], *q_s, t,*this);  
}
//added 11.04.2015
CCStateSpace::CCTurn CCStateSpace::CC_SN(CCStateSpace::State *q_s, double t)
{
    return CCStateSpace::CCTurn(ccTurnType[1], *q_s, t,*this);
}
//--------------------------------End of Primitive CCTurn-------------------------------------
//--------------------------------------------------------------------------------------------
//added 11.04.2015, verified
//CSC,8
//1
bool CCStateSpace::LPSPLP(CCStateSpace::State *q_s,CCStateSpace::State *q_g, vector<CCStateSpace::CCTurn> *vec_cct, double *t, double *u, double *v, double *err, CCStateSpace::CCPathType *type)
{
    double uu = dist2d(q_g->y_o_ln - q_s->y_o_lp, q_g->x_o_ln - q_s->x_o_lp)-2 * this->r*sin(this->miu);
    if (uu > -CC_EPS)
    {
        *u=uu;
        *type=LpSpLp;
        double gamma = mod2pi2(atan2(q_g->y_o_ln - q_s->y_o_lp, q_g->x_o_ln - q_s->x_o_lp));
        *t = mod2pi2(gamma-q_s->theta);
        *v = mod2pi2(q_g->theta-gamma);

        if(!vec_cct->empty()) vec_cct->clear();
        vec_cct->push_back(this->CC_LP(q_s, *t));
        if(*u>CC_EPS) vec_cct->push_back(this->CC_SP(&(vec_cct->back()).q_g, *u));
        vec_cct->push_back(this->CC_LP(&(vec_cct->back()).q_g, *v));

        *err=this->distance(&(vec_cct->back()).q_g, q_g);
        return true;
    }
    else
        return false;
}
//2
bool CCStateSpace::LNSNLN(CCStateSpace::State *q_s,CCStateSpace::State *q_g, vector<CCStateSpace::CCTurn> *vec_cct, double *t, double *u, double *v, double *err, CCStateSpace::CCPathType *type)
{
    double uu = dist2d(q_g->y_o_lp - q_s->y_o_ln, q_g->x_o_lp - q_s->x_o_ln)-2 * this->r*sin(this->miu);
    if (uu > -CC_EPS)
    {
        *u=uu;
        *type=LnSnLn;
        double gamma = mod2pi2(atan2(q_g->y_o_lp - q_s->y_o_ln, q_g->x_o_lp - q_s->x_o_ln));
        *t = abs(mod2pi3(pi+gamma-q_s->theta));
        *v = abs(mod2pi3(q_g->theta-gamma-pi));

        if(!vec_cct->empty()) vec_cct->clear();
        vec_cct->push_back(this->CC_LN(q_s, *t));
        if(*u>CC_EPS) vec_cct->push_back(this->CC_SN(&(vec_cct->back()).q_g, *u));
        vec_cct->push_back(this->CC_LN(&(vec_cct->back()).q_g, *v));

        *err=this->distance(&(vec_cct->back()).q_g, q_g);
        return true;
    }
    else
        return false;
}
//3
bool CCStateSpace::RPSPRP(CCStateSpace::State *q_s,CCStateSpace::State *q_g, vector<CCStateSpace::CCTurn> *vec_cct, double *t, double *u, double *v, double *err, CCStateSpace::CCPathType *type)
{
    double uu = dist2d(q_g->y_o_rn - q_s->y_o_rp, q_g->x_o_rn - q_s->x_o_rp)-2 * this->r*sin(this->miu);
    if (uu > -CC_EPS)
    {
        *u=uu;
        *type=RpSpRp;
        double gamma = mod2pi2(atan2(q_g->y_o_rn - q_s->y_o_rp, q_g->x_o_rn - q_s->x_o_rp));
        *t = abs(mod2pi3(gamma-q_s->theta));
        *v = abs(mod2pi3(q_g->theta-gamma));

        if(!vec_cct->empty()) vec_cct->clear();
        vec_cct->push_back(this->CC_RP(q_s, *t));
        if(*u>CC_EPS) vec_cct->push_back(this->CC_SP(&(vec_cct->back()).q_g, *u));
        vec_cct->push_back(this->CC_RP(&(vec_cct->back()).q_g, *v));

        *err=this->distance(&(vec_cct->back()).q_g, q_g);
        return true;
    }
    else
        return false;
}
//4
bool CCStateSpace::RNSNRN(CCStateSpace::State *q_s,CCStateSpace::State *q_g, vector<CCStateSpace::CCTurn> *vec_cct, double *t, double *u, double *v, double *err, CCStateSpace::CCPathType *type)
{
    double uu = dist2d(q_g->y_o_rp - q_s->y_o_rn, q_g->x_o_rp - q_s->x_o_rn)-2 * this->r*sin(this->miu);
    if (uu > -CC_EPS)
    {
        *u=uu;
        *type=RnSnRn;
        double gamma = mod2pi2(atan2(q_g->y_o_rp - q_s->y_o_rn, q_g->x_o_rp - q_s->x_o_rn));
        *t = mod2pi2(pi+gamma-q_s->theta);
        *v = mod2pi2(q_g->theta-gamma-pi);

        if(!vec_cct->empty()) vec_cct->clear();
        vec_cct->push_back(this->CC_RN(q_s, *t));
        if(*u>CC_EPS) vec_cct->push_back(this->CC_SN(&(vec_cct->back()).q_g, *u));
        vec_cct->push_back(this->CC_RN(&(vec_cct->back()).q_g, *v));

        *err=this->distance(&(vec_cct->back()).q_g, q_g);
        return true;
    }
    else
        return false;
}
//5
bool CCStateSpace::LPSPRP(CCStateSpace::State *q_s,CCStateSpace::State *q_g, vector<CCStateSpace::CCTurn> *vec_cct, double *t, double *u, double *v, double *err, CCStateSpace::CCPathType *type)
{
    double uu = sqrt(dist2(q_g->y_o_rn - q_s->y_o_lp, q_g->x_o_rn - q_s->x_o_lp)-4*this->r*this->r*cos(this->miu)*cos(this->miu))-2 * this->r*sin(this->miu);
    if (uu > -CC_EPS)
    {
        *u=uu;
        *type=LpSpRp;
        double gamma = mod2pi2(atan2(q_g->y_o_rn - q_s->y_o_lp, q_g->x_o_rn - q_s->x_o_lp)+asin(2*this->r*sin(this->miu)/dist2d(q_g->y_o_rn - q_s->y_o_lp, q_g->x_o_rn - q_s->x_o_lp)));
        *t = mod2pi2(gamma-q_s->theta);
        *v = abs(mod2pi3(q_g->theta-gamma));

        if(!vec_cct->empty()) vec_cct->clear();
        vec_cct->push_back(this->CC_LP(q_s, *t));
        if(*u>CC_EPS) vec_cct->push_back(this->CC_SP(&(vec_cct->back()).q_g, *u));
        vec_cct->push_back(this->CC_RP(&(vec_cct->back()).q_g, *v));

        *err=this->distance(&(vec_cct->back()).q_g, q_g);
        return true;
    }
    else
        return false;
}
//6
bool CCStateSpace::LNSNRN(CCStateSpace::State *q_s,CCStateSpace::State *q_g, vector<CCStateSpace::CCTurn> *vec_cct, double *t, double *u, double *v, double *err, CCStateSpace::CCPathType *type)
{
    double uu = sqrt(dist2(q_g->y_o_rp - q_s->y_o_ln, q_g->x_o_rp - q_s->x_o_ln)-4*this->r*this->r*cos(this->miu)*cos(this->miu))-2 * this->r*sin(this->miu);
    if (uu > -CC_EPS)
    {
        *u=uu;
        *type=LnSnRn;
        double gamma = mod2pi2(atan2(q_g->y_o_rp - q_s->y_o_ln, q_g->x_o_rp - q_s->x_o_ln) - asin(2*this->r*sin(this->miu)/dist2d(q_g->y_o_rp - q_s->y_o_ln, q_g->x_o_rp - q_s->x_o_ln)));
        *t = abs(mod2pi3(gamma+pi-q_s->theta));
        *v = mod2pi2(q_g->theta-gamma-pi);

        if(!vec_cct->empty()) vec_cct->clear();
        vec_cct->push_back(this->CC_LN(q_s, *t));
        if(*u>CC_EPS) vec_cct->push_back(this->CC_SN(&(vec_cct->back()).q_g, *u));
        vec_cct->push_back(this->CC_RN(&(vec_cct->back()).q_g, *v));

        *err=this->distance(&(vec_cct->back()).q_g, q_g);
        return true;
    }
    else
        return false;
}
//7
bool CCStateSpace::RPSPLP(CCStateSpace::State *q_s,CCStateSpace::State *q_g, vector<CCStateSpace::CCTurn> *vec_cct, double *t, double *u, double *v, double *err, CCStateSpace::CCPathType *type)
{
    double uu = sqrt(dist2(q_g->y_o_ln - q_s->y_o_rp, q_g->x_o_ln - q_s->x_o_rp)-4*this->r*this->r*cos(this->miu)*cos(this->miu))-2 * this->r*sin(this->miu);
    if (uu > -CC_EPS)
    {
        *u=uu;
        *type=RpSpLp;
        double gamma = mod2pi2(atan2(q_g->y_o_ln - q_s->y_o_rp, q_g->x_o_ln - q_s->x_o_rp) - asin(2*this->r*sin(this->miu)/dist2d(q_g->y_o_ln - q_s->y_o_rp, q_g->x_o_ln - q_s->x_o_rp)));
        *t = abs(mod2pi3(gamma-q_s->theta));
        *v = mod2pi2(q_g->theta-gamma);

        if(!vec_cct->empty()) vec_cct->clear();
        vec_cct->push_back(this->CC_RP(q_s, *t));
        if(*u>CC_EPS) vec_cct->push_back(this->CC_SP(&(vec_cct->back()).q_g, *u));
        vec_cct->push_back(this->CC_LP(&(vec_cct->back()).q_g, *v));

        *err=this->distance(&(vec_cct->back()).q_g, q_g);
        return true;
    }
    else
        return false;
}
//8
bool CCStateSpace::RNSNLN(CCStateSpace::State *q_s,CCStateSpace::State *q_g, vector<CCStateSpace::CCTurn> *vec_cct, double *t, double *u, double *v, double *err, CCStateSpace::CCPathType *type)
{
    double uu = sqrt(dist2(q_g->y_o_lp - q_s->y_o_rn, q_g->x_o_lp - q_s->x_o_rn)-4*this->r*this->r*cos(this->miu)*cos(this->miu))-2 * this->r*sin(this->miu);
    if (uu > -CC_EPS)
    {
        *u=uu;
        *type=RnSnLn;
        double gamma = mod2pi2(atan2(q_g->y_o_lp - q_s->y_o_rn, q_g->x_o_lp - q_s->x_o_rn) + asin(2*this->r*sin(this->miu)/dist2d(q_g->y_o_lp - q_s->y_o_rn, q_g->x_o_lp - q_s->x_o_rn)));
        *t = mod2pi2(gamma+pi-q_s->theta);
        *v = abs(mod2pi3(q_g->theta-gamma-pi));

        if(!vec_cct->empty()) vec_cct->clear();
        vec_cct->push_back(this->CC_RN(q_s, *t));
        if(*u>CC_EPS) vec_cct->push_back(this->CC_SN(&(vec_cct->back()).q_g, *u));
        vec_cct->push_back(this->CC_LN(&(vec_cct->back()).q_g, *v));

        *err=this->distance(&(vec_cct->back()).q_g, q_g);
        return true;
    }
    else
        return false;
}

//CCC,12-------------------------------------------------------
//9
bool CCStateSpace::LPRNLP(CCStateSpace::State *q_s,CCStateSpace::State *q_g, vector<CCStateSpace::CCTurn> *vec_cct, double *t, double *u, double *v, double *err, CCStateSpace::CCPathType *type)
{
    double l_max = 4*this->r*cos(this->miu);
    double l= dist2d(q_g->y_o_ln - q_s->y_o_lp, q_g->x_o_ln - q_s->x_o_lp);
    if (l <= l_max + CC_EPS)
    {
        *type=LpRnLp;
        double gamma=mod2pi2(atan2(q_g->y_o_ln - q_s->y_o_lp, q_g->x_o_ln - q_s->x_o_lp));
        double beta=(l>l_max-CC_EPS)?0:acos(l/l_max);
        *t=mod2pi2(gamma+beta+half_pi-q_s->theta);
        *u=mod2pi2(pi-2*beta);
        *v=mod2pi2(q_g->theta-gamma+beta+half_pi);

        if(!vec_cct->empty()) vec_cct->clear();
        vec_cct->push_back(this->CC_LP(q_s, *t));
        vec_cct->push_back(this->CC_RN(&(vec_cct->back()).q_g, *u));
        vec_cct->push_back(this->CC_LP(&(vec_cct->back()).q_g, *v));

        *err=this->distance(&(vec_cct->back()).q_g, q_g);
        return true;
    }
    else
        return false;
}
//10
bool CCStateSpace::LNRPLN(CCStateSpace::State *q_s,CCStateSpace::State *q_g, vector<CCStateSpace::CCTurn> *vec_cct, double *t, double *u, double *v, double *err, CCStateSpace::CCPathType *type)
{
    double l_max = 4*this->r*cos(this->miu);
    double l= dist2d(q_g->y_o_lp - q_s->y_o_ln, q_g->x_o_lp - q_s->x_o_ln);
    if (l <= l_max + CC_EPS)
    {
        *type=LnRpLn;
        double gamma=mod2pi2(atan2(q_g->y_o_lp - q_s->y_o_ln, q_g->x_o_lp - q_s->x_o_ln));
        double beta=(l>l_max-CC_EPS)?0:acos(l/l_max);
        *t=abs(mod2pi3(gamma-beta+half_pi-q_s->theta));
        *u=abs(mod2pi3(2*beta-pi));
        *v=abs(mod2pi3(q_g->theta-gamma-beta+half_pi));

        if(!vec_cct->empty()) vec_cct->clear();
        vec_cct->push_back(this->CC_LN(q_s, *t));
        vec_cct->push_back(this->CC_RP(&(vec_cct->back()).q_g, *u));
        vec_cct->push_back(this->CC_LN(&(vec_cct->back()).q_g, *v));

        *err=this->distance(&(vec_cct->back()).q_g, q_g);
        return true;
    }
    else
        return false;
}
//11
bool CCStateSpace::RPLNRP(CCStateSpace::State *q_s,CCStateSpace::State *q_g, vector<CCStateSpace::CCTurn> *vec_cct, double *t, double *u, double *v, double *err, CCStateSpace::CCPathType *type)
{
    double l_max = 4*this->r*cos(this->miu);
    double l= dist2d(q_g->y_o_rn - q_s->y_o_rp, q_g->x_o_rn - q_s->x_o_rp);
    if (l <= l_max + CC_EPS)
    {
        *type=RpLnRp;
        double gamma=mod2pi2(atan2(q_g->y_o_rn - q_s->y_o_rp, q_g->x_o_rn - q_s->x_o_rp));
        double beta=(l>l_max-CC_EPS)?0:acos(l/l_max);
        *t=abs(mod2pi3(gamma-beta-half_pi-q_s->theta));
        *u=abs(mod2pi3(2*beta+pi));
        *v=abs(mod2pi3(q_g->theta-gamma-beta-half_pi));

        if(!vec_cct->empty()) vec_cct->clear();
        vec_cct->push_back(this->CC_RP(q_s, *t));
        vec_cct->push_back(this->CC_LN(&(vec_cct->back()).q_g, *u));
        vec_cct->push_back(this->CC_RP(&(vec_cct->back()).q_g, *v));

        *err=this->distance(&(vec_cct->back()).q_g, q_g);
        return true;
    }
    else
        return false;
}
//12
bool CCStateSpace::RNLPRN(CCStateSpace::State *q_s,CCStateSpace::State *q_g, vector<CCStateSpace::CCTurn> *vec_cct, double *t, double *u, double *v, double *err, CCStateSpace::CCPathType *type)
{
    double l_max = 4*this->r*cos(this->miu);
    double l= dist2d(q_g->y_o_rp - q_s->y_o_rn, q_g->x_o_rp - q_s->x_o_rn);
    if (l <= l_max + CC_EPS)
    {
        *type=RnLpRn;
        double gamma=mod2pi2(atan2(q_g->y_o_rp - q_s->y_o_rn, q_g->x_o_rp - q_s->x_o_rn));
        double beta=(l>l_max-CC_EPS)?0:acos(l/l_max);
        *t=mod2pi2(gamma+beta-half_pi-q_s->theta);
        *u=mod2pi2(pi-2*beta);
        *v=mod2pi2(q_g->theta-gamma+beta-half_pi);

        if(!vec_cct->empty()) vec_cct->clear();
        vec_cct->push_back(this->CC_RN(q_s, *t));
        vec_cct->push_back(this->CC_LP(&(vec_cct->back()).q_g, *u));
        vec_cct->push_back(this->CC_RN(&(vec_cct->back()).q_g, *v));

        *err=this->distance(&(vec_cct->back()).q_g, q_g);
        return true;
    }
    else
        return false;
}
//13
bool CCStateSpace::LPRNLN(CCStateSpace::State *q_s,CCStateSpace::State *q_g, vector<CCStateSpace::CCTurn> *vec_cct, double *t, double *u, double *v, double *err, CCStateSpace::CCPathType *type)
{
    double l_max = 2*this->r*(1+cos(this->miu));
    double l= dist2d(q_g->y_o_lp - q_s->y_o_lp, q_g->x_o_lp - q_s->x_o_lp);
    if (l <= l_max + CC_EPS)
    {
        *type=LpRnLn;
        double gamma=mod2pi2(atan2(q_g->y_o_lp - q_s->y_o_lp, q_g->x_o_lp - q_s->x_o_lp));
        double l1=l/2-2*this->r*this->r*sin(this->miu)*sin(this->miu)/l;
        double l2=l-l1;
        double beta1=(l>l_max-CC_EPS)?0:acos(l1/(2*this->r*cos(this->miu)));
        double beta2=(l>l_max-CC_EPS)?0:acos(l2/(2*this->r));

        *t=mod2pi2(gamma+beta1+half_pi-q_s->theta);
        *u=mod2pi2(pi-beta1-beta2-this->miu);
        *v=abs(mod2pi3(q_g->theta-gamma+beta2+this->miu+half_pi));

        if(!vec_cct->empty()) vec_cct->clear();
        vec_cct->push_back(this->CC_LP(q_s, *t));
        vec_cct->push_back(this->CC_RN(&(vec_cct->back()).q_g, *u));
        vec_cct->push_back(this->CC_LN(&(vec_cct->back()).q_g, *v));

        *err=this->distance(&(vec_cct->back()).q_g, q_g);
        return true;
    }
    else
        return false;
}
//14
bool CCStateSpace::LNRPLP(CCStateSpace::State *q_s,CCStateSpace::State *q_g, vector<CCStateSpace::CCTurn> *vec_cct, double *t, double *u, double *v, double *err, CCStateSpace::CCPathType *type)
{
    double l_max = 2*this->r*(1+cos(this->miu));
    double l= dist2d(q_g->y_o_ln - q_s->y_o_ln, q_g->x_o_ln - q_s->x_o_ln);
    if (l <= l_max + CC_EPS)
    {
        *type=LnRpLp;
        double gamma=mod2pi2(atan2(q_g->y_o_ln - q_s->y_o_ln, q_g->x_o_ln - q_s->x_o_ln));
        double l1=l/2-2*this->r*this->r*sin(this->miu)*sin(this->miu)/l;
        double l2=l-l1;
        double beta1=(l>l_max-CC_EPS)?0:acos(l1/(2*this->r*cos(this->miu)));
        double beta2=(l>l_max-CC_EPS)?0:acos(l2/(2*this->r));

        *t=abs(mod2pi3(gamma-beta1+half_pi-q_s->theta));
        *u=abs(mod2pi3(this->miu+beta1+beta2+pi));
        *v=mod2pi2(q_g->theta-gamma-beta2-this->miu+half_pi);

        if(!vec_cct->empty()) vec_cct->clear();
        vec_cct->push_back(this->CC_LN(q_s, *t));
        vec_cct->push_back(this->CC_RP(&(vec_cct->back()).q_g, *u));
        vec_cct->push_back(this->CC_LP(&(vec_cct->back()).q_g, *v));

        *err=this->distance(&(vec_cct->back()).q_g, q_g);
        return true;
    }
    else
        return false;
}
//15
bool CCStateSpace::RPLNRN(CCStateSpace::State *q_s,CCStateSpace::State *q_g, vector<CCStateSpace::CCTurn> *vec_cct, double *t, double *u, double *v, double *err, CCStateSpace::CCPathType *type)
{
    double l_max = 2*this->r*(1+cos(this->miu));
    double l= dist2d(q_g->y_o_rp - q_s->y_o_rp, q_g->x_o_rp - q_s->x_o_rp);
    if (l <= l_max + CC_EPS)
    {
        *type=RpLnRn;
        double gamma=mod2pi2(atan2(q_g->y_o_rp - q_s->y_o_rp, q_g->x_o_rp - q_s->x_o_rp));
        double l1=l/2-2*this->r*this->r*sin(this->miu)*sin(this->miu)/l;
        double l2=l-l1;
        double beta1=(l>l_max-CC_EPS)?0:acos(l1/(2*this->r*cos(this->miu)));
        double beta2=(l>l_max-CC_EPS)?0:acos(l2/(2*this->r));

        *t=abs(mod2pi3(gamma-beta1-half_pi-q_s->theta));
        *u=abs(mod2pi3(this->miu+beta1+beta2+pi));
        *v=mod2pi2(q_g->theta-gamma-beta2-this->miu-half_pi);

        if(!vec_cct->empty()) vec_cct->clear();
        vec_cct->push_back(this->CC_RP(q_s, *t));
        vec_cct->push_back(this->CC_LN(&(vec_cct->back()).q_g, *u));
        vec_cct->push_back(this->CC_RN(&(vec_cct->back()).q_g, *v));

        *err=this->distance(&(vec_cct->back()).q_g, q_g);
        return true;
    }
    else
        return false;
}
//16
bool CCStateSpace::RNLPRP(CCStateSpace::State *q_s,CCStateSpace::State *q_g, vector<CCStateSpace::CCTurn> *vec_cct, double *t, double *u, double *v, double *err, CCStateSpace::CCPathType *type)
{
    double l_max = 2*this->r*(1+cos(this->miu));
    double l= dist2d(q_g->y_o_rn - q_s->y_o_rn, q_g->x_o_rn - q_s->x_o_rn);
    if (l <= l_max + CC_EPS)
    {
        *type=RnLpRp;
        double gamma=mod2pi2(atan2(q_g->y_o_rn - q_s->y_o_rn, q_g->x_o_rn - q_s->x_o_rn));
        double l1=l/2-2*this->r*this->r*sin(this->miu)*sin(this->miu)/l;
        double l2=l-l1;
        double beta1=(l>l_max-CC_EPS)?0:acos(l1/(2*this->r*cos(this->miu)));
        double beta2=(l>l_max-CC_EPS)?0:acos(l2/(2*this->r));

        *t=mod2pi2(gamma+beta1-half_pi-q_s->theta);
        *u=mod2pi2(pi-this->miu-beta1-beta2);
        *v=abs(mod2pi3(q_g->theta-gamma+beta2+this->miu-half_pi));

        if(!vec_cct->empty()) vec_cct->clear();
        vec_cct->push_back(this->CC_RN(q_s, *t));
        vec_cct->push_back(this->CC_LP(&(vec_cct->back()).q_g, *u));
        vec_cct->push_back(this->CC_RP(&(vec_cct->back()).q_g, *v));

        *err=this->distance(&(vec_cct->back()).q_g, q_g);
        return true;
    }
    else
        return false;
}
//17
bool CCStateSpace::LNRNLP(CCStateSpace::State *q_s,CCStateSpace::State *q_g, vector<CCStateSpace::CCTurn> *vec_cct, double *t, double *u, double *v, double *err, CCStateSpace::CCPathType *type)
{
    double l_max = 2*this->r*(1+cos(this->miu));
    double l= dist2d(q_g->y_o_ln - q_s->y_o_ln, q_g->x_o_ln - q_s->x_o_ln);
    if (l <= l_max + CC_EPS)
    {
        *type=LnRnLp;
        double gamma=mod2pi2(atan2(q_g->y_o_ln - q_s->y_o_ln, q_g->x_o_ln - q_s->x_o_ln));
        double l1=l/2+2*this->r*this->r*sin(this->miu)*sin(this->miu)/l;
        double l2=l-l1;
        double beta1=(l>l_max-CC_EPS)?0:acos(l2/(2*this->r));
        double beta2=(l>l_max-CC_EPS)?0:acos(l1/(2*this->r*cos(this->miu)));

        *t=abs(mod2pi3(gamma-beta1+this->miu+half_pi-q_s->theta));
        *u=mod2pi2(beta1+beta2-pi-this->miu);
        *v=mod2pi2(q_g->theta-gamma-beta2+half_pi);

        if(!vec_cct->empty()) vec_cct->clear();
        vec_cct->push_back(this->CC_LN(q_s, *t));
        vec_cct->push_back(this->CC_RN(&(vec_cct->back()).q_g, *u));
        vec_cct->push_back(this->CC_LP(&(vec_cct->back()).q_g, *v));

        *err=this->distance(&(vec_cct->back()).q_g, q_g);
        return true;
    }
    else
        return false;
}
//18
bool CCStateSpace::LPRPLN(CCStateSpace::State *q_s,CCStateSpace::State *q_g, vector<CCStateSpace::CCTurn> *vec_cct, double *t, double *u, double *v, double *err, CCStateSpace::CCPathType *type)
{
    double l_max = 2*this->r*(1+cos(this->miu));
    double l= dist2d(q_g->y_o_lp - q_s->y_o_lp, q_g->x_o_lp - q_s->x_o_lp);
    if (l <= l_max + CC_EPS)
    {
        *type=LpRpLn;
        double gamma=mod2pi2(atan2(q_g->y_o_lp - q_s->y_o_lp, q_g->x_o_lp - q_s->x_o_lp));
        double l1=l/2+2*this->r*this->r*sin(this->miu)*sin(this->miu)/l;
        double l2=l-l1;
        double beta1=(l>l_max-CC_EPS)?0:acos(l2/(2*this->r));
        double beta2=(l>l_max-CC_EPS)?0:acos(l1/(2*this->r*cos(this->miu)));

        *t=mod2pi2(gamma-beta1-this->miu+half_pi-q_s->theta);
        *u=abs(mod2pi3(beta1+beta2-pi+this->miu));
        *v=abs(mod2pi3(q_g->theta-gamma-beta2+half_pi));

        if(!vec_cct->empty()) vec_cct->clear();
        vec_cct->push_back(this->CC_LP(q_s, *t));
        vec_cct->push_back(this->CC_RP(&(vec_cct->back()).q_g, *u));
        vec_cct->push_back(this->CC_LN(&(vec_cct->back()).q_g, *v));

        *err=this->distance(&(vec_cct->back()).q_g, q_g);
        return true;
    }
    else
        return false;
}
//19
bool CCStateSpace::RNLNRP(CCStateSpace::State *q_s,CCStateSpace::State *q_g, vector<CCStateSpace::CCTurn> *vec_cct, double *t, double *u, double *v, double *err, CCStateSpace::CCPathType *type)
{
    double l_max = 2*this->r*(1+cos(this->miu));
    double l= dist2d(q_g->y_o_rn - q_s->y_o_rn, q_g->x_o_rn - q_s->x_o_rn);
    if (l <= l_max + CC_EPS)
    {
        *type=RnLnRp;
        double gamma=mod2pi2(atan2(q_g->y_o_rn - q_s->y_o_rn, q_g->x_o_rn - q_s->x_o_rn));
        double l1=l/2+2*this->r*this->r*sin(this->miu)*sin(this->miu)/l;
        double l2=l-l1;
        double beta1=(l>l_max-CC_EPS)?0:acos(l2/(2*this->r));
        double beta2=(l>l_max-CC_EPS)?0:acos(l1/(2*this->r*cos(this->miu)));

        *t=mod2pi2(gamma-beta1-this->miu-half_pi-q_s->theta);
        *u=abs(mod2pi3(beta1+beta2+pi+this->miu));
        *v=abs(mod2pi3(q_g->theta-gamma-beta2-half_pi));

        if(!vec_cct->empty()) vec_cct->clear();
        vec_cct->push_back(this->CC_RN(q_s, *t));
        vec_cct->push_back(this->CC_LN(&(vec_cct->back()).q_g, *u));
        vec_cct->push_back(this->CC_RP(&(vec_cct->back()).q_g, *v));

        *err=this->distance(&(vec_cct->back()).q_g, q_g);
        return true;
    }
    else
        return false;
}
//20
bool CCStateSpace::RPLPRN(CCStateSpace::State *q_s,CCStateSpace::State *q_g, vector<CCStateSpace::CCTurn> *vec_cct, double *t, double *u, double *v, double *err, CCStateSpace::CCPathType *type)
{
    double l_max = 2*this->r*(1+cos(this->miu));
    double l= dist2d(q_g->y_o_rp - q_s->y_o_rp, q_g->x_o_rp - q_s->x_o_rp);
    if (l <= l_max + CC_EPS)
    {
        *type=RpLpRn;
        double gamma=mod2pi2(atan2(q_g->y_o_rp - q_s->y_o_rp, q_g->x_o_rp - q_s->x_o_rp));
        double l1=l/2+2*this->r*this->r*sin(this->miu)*sin(this->miu)/l;
        double l2=l-l1;
        double beta1=(l>l_max-CC_EPS)?0:acos(l2/(2*this->r));
        double beta2=(l>l_max-CC_EPS)?0:acos(l1/(2*this->r*cos(this->miu)));

        *t=abs(mod2pi3(gamma+beta1+this->miu-half_pi-q_s->theta));
        *u=mod2pi2(pi-beta1-beta2-this->miu);
        *v=mod2pi2(q_g->theta-gamma+beta2-half_pi);

        if(!vec_cct->empty()) vec_cct->clear();
        vec_cct->push_back(this->CC_RP(q_s, *t));
        vec_cct->push_back(this->CC_LP(&(vec_cct->back()).q_g, *u));
        vec_cct->push_back(this->CC_RN(&(vec_cct->back()).q_g, *v));

        *err=this->distance(&(vec_cct->back()).q_g, q_g);
        return true;
    }
    else
        return false;
}
//CCCC,8-------------------------------------------------------
//21
bool CCStateSpace::LPRPLNRN(CCStateSpace::State *q_s,CCStateSpace::State *q_g, vector<CCStateSpace::CCTurn> *vec_cct, double *t, double *u, double *v, double *err, CCStateSpace::CCPathType *type)
{
    double l_max = 2*(2+cos(this->miu))*this->r;
    double l= dist2d(q_g->y_o_rp - q_s->y_o_lp, q_g->x_o_rp - q_s->x_o_lp);
    if (l <= l_max + CC_EPS)
    {
        *type=LpRpLnRn;
        double gamma=mod2pi2(atan2(q_g->y_o_rp - q_s->y_o_lp, q_g->x_o_rp - q_s->x_o_lp));
        double beta= (l>l_max-CC_EPS)?0:acos((l-2*this->r*cos(this->miu))/(4*this->r));

        *t=mod2pi2(gamma-beta-this->miu+half_pi-q_s->theta);
        *u=abs(mod2pi3(beta+this->miu+pi));
        *v=mod2pi2(q_g->theta-gamma-beta-this->miu-half_pi);

        if(!vec_cct->empty()) vec_cct->clear();
        vec_cct->push_back(this->CC_LP(q_s, *t));
        vec_cct->push_back(this->CC_RP(&(vec_cct->back()).q_g, *u));
        vec_cct->push_back(this->CC_LN(&(vec_cct->back()).q_g, *u));
        vec_cct->push_back(this->CC_RN(&(vec_cct->back()).q_g, *v));

        *err=this->distance(&(vec_cct->back()).q_g, q_g);
        return true;
    }
    else
        return false;
}
//22
bool CCStateSpace::LNRNLPRP(CCStateSpace::State *q_s,CCStateSpace::State *q_g, vector<CCStateSpace::CCTurn> *vec_cct, double *t, double *u, double *v, double *err, CCStateSpace::CCPathType *type)
{
    double l_max = 2*(2+cos(this->miu))*this->r;
    double l= dist2d(q_g->y_o_rn - q_s->y_o_ln, q_g->x_o_rn - q_s->x_o_ln);
    if (l <= l_max + CC_EPS)
    {
        *type=LnRnLpRp;
        double gamma=mod2pi2(atan2(q_g->y_o_rn - q_s->y_o_ln, q_g->x_o_rn - q_s->x_o_ln));
        double beta= (l>l_max-CC_EPS)?0:acos((l-2*this->r*cos(this->miu))/(4*this->r));

        *t=abs(mod2pi3(gamma+beta+this->miu+half_pi-q_s->theta));
        *u=mod2pi2(pi-beta-this->miu);
        *v=abs(mod2pi3(q_g->theta-gamma+beta+this->miu-half_pi));

        if(!vec_cct->empty()) vec_cct->clear();
        vec_cct->push_back(this->CC_LN(q_s, *t));
        vec_cct->push_back(this->CC_RN(&(vec_cct->back()).q_g, *u));
        vec_cct->push_back(this->CC_LP(&(vec_cct->back()).q_g, *u));
        vec_cct->push_back(this->CC_RP(&(vec_cct->back()).q_g, *v));

        *err=this->distance(&(vec_cct->back()).q_g, q_g);
        return true;
    }
    else
        return false;
}
//23
bool CCStateSpace::RPLPRNLN(CCStateSpace::State *q_s,CCStateSpace::State *q_g, vector<CCStateSpace::CCTurn> *vec_cct, double *t, double *u, double *v, double *err, CCStateSpace::CCPathType *type)
{
    double l_max = 2*(2+cos(this->miu))*this->r;
    double l= dist2d(q_g->y_o_lp - q_s->y_o_rp, q_g->x_o_lp - q_s->x_o_rp);
    if (l <= l_max + CC_EPS)
    {
        *type=RpLpRnLn;
        double gamma=mod2pi2(atan2(q_g->y_o_lp - q_s->y_o_rp, q_g->x_o_lp - q_s->x_o_rp));
        double beta= (l>l_max-CC_EPS)?0:acos((l-2*this->r*cos(this->miu))/(4*this->r));

        *t=abs(mod2pi3(gamma+beta+this->miu-half_pi-q_s->theta));
        *u=mod2pi2(pi-beta-this->miu);
        *v=abs(mod2pi3(q_g->theta-gamma+beta+this->miu+half_pi));

        if(!vec_cct->empty()) vec_cct->clear();
        vec_cct->push_back(this->CC_RP(q_s, *t));
        vec_cct->push_back(this->CC_LP(&(vec_cct->back()).q_g, *u));
        vec_cct->push_back(this->CC_RN(&(vec_cct->back()).q_g, *u));
        vec_cct->push_back(this->CC_LN(&(vec_cct->back()).q_g, *v));

        *err=this->distance(&(vec_cct->back()).q_g, q_g);
        return true;
    }
    else
        return false;
}
//24
bool CCStateSpace::RNLNRPLP(CCStateSpace::State *q_s,CCStateSpace::State *q_g, vector<CCStateSpace::CCTurn> *vec_cct, double *t, double *u, double *v, double *err, CCStateSpace::CCPathType *type)
{
    double l_max = 2*(2+cos(this->miu))*this->r;
    double l= dist2d(q_g->y_o_ln - q_s->y_o_rn, q_g->x_o_ln - q_s->x_o_rn);
    if (l <= l_max + CC_EPS)
    {
        *type=RnLnRpLp;
        double gamma=mod2pi2(atan2(q_g->y_o_ln - q_s->y_o_rn, q_g->x_o_ln - q_s->x_o_rn));
        double beta= (l>l_max-CC_EPS)?0:acos((l-2*this->r*cos(this->miu))/(4*this->r));

        *t=mod2pi2(gamma-beta-this->miu-half_pi-q_s->theta);
        *u=abs(mod2pi3(pi+beta+this->miu));
        *v=mod2pi2(q_g->theta-gamma-beta-this->miu+half_pi);

        if(!vec_cct->empty()) vec_cct->clear();
        vec_cct->push_back(this->CC_RN(q_s, *t));
        vec_cct->push_back(this->CC_LN(&(vec_cct->back()).q_g, *u));
        vec_cct->push_back(this->CC_RP(&(vec_cct->back()).q_g, *u));
        vec_cct->push_back(this->CC_LP(&(vec_cct->back()).q_g, *v));

        *err=this->distance(&(vec_cct->back()).q_g, q_g);
        return true;
    }
    else
        return false;
}
//25
bool CCStateSpace::LPRNLNRP(CCStateSpace::State *q_s,CCStateSpace::State *q_g, vector<CCStateSpace::CCTurn> *vec_cct, double *t, double *u, double *v, double *err, CCStateSpace::CCPathType *type)
{
    double l_max = 2*(1+2*cos(this->miu))*this->r;
    double l= dist2d(q_g->y_o_rn - q_s->y_o_lp, q_g->x_o_rn - q_s->x_o_lp);
    if (l <= l_max + CC_EPS)
    {
        *type=LpRnLnRp;
        double gamma=mod2pi2(atan2(q_g->y_o_rn - q_s->y_o_lp, q_g->x_o_rn - q_s->x_o_lp));
        double l1=l/4+this->r*this->r*(4*cos(this->miu)*cos(this->miu)-1)/l;
        double l2=l/2-l1;
        double beta1=(l>l_max-CC_EPS)?0:acos(l1/(2*this->r*cos(this->miu)));
        double beta2=(l>l_max-CC_EPS)?0:acos(l2/(this->r));

        *t=mod2pi2(gamma+beta1+half_pi-q_s->theta);
        *u=mod2pi2(pi-beta1-beta2-this->miu);
        *v=abs(mod2pi3(q_g->theta-gamma-beta1-half_pi));

        if(!vec_cct->empty()) vec_cct->clear();
        vec_cct->push_back(this->CC_LP(q_s, *t));
        vec_cct->push_back(this->CC_RN(&(vec_cct->back()).q_g, *u));
        vec_cct->push_back(this->CC_LN(&(vec_cct->back()).q_g, *u));
        vec_cct->push_back(this->CC_RP(&(vec_cct->back()).q_g, *v));

        *err=this->distance(&(vec_cct->back()).q_g, q_g);
        return true;
    }
    else
        return false;
}
//26
bool CCStateSpace::LNRPLPRN(CCStateSpace::State *q_s,CCStateSpace::State *q_g, vector<CCStateSpace::CCTurn> *vec_cct, double *t, double *u, double *v, double *err, CCStateSpace::CCPathType *type)
{
    double l_max = 2*(1+2*cos(this->miu))*this->r;
    double l= dist2d(q_g->y_o_rp - q_s->y_o_ln, q_g->x_o_rp - q_s->x_o_ln);
    if (l <= l_max + CC_EPS)
    {
        *type=LnRpLpRn;
        double gamma=mod2pi2(atan2(q_g->y_o_rp - q_s->y_o_ln, q_g->x_o_rp - q_s->x_o_ln));
        double l1=l/4+this->r*this->r*(4*cos(this->miu)*cos(this->miu)-1)/l;
        double l2=l/2-l1;
        double beta1=(l>l_max-CC_EPS)?0:acos(l1/(2*this->r*cos(this->miu)));
        double beta2=(l>l_max-CC_EPS)?0:acos(l2/(this->r));

        *t=abs(mod2pi3(gamma-beta1+half_pi-q_s->theta));
        *u=mod2pi2(pi-beta1-beta2-this->miu);
        *v=mod2pi2(q_g->theta-gamma+beta1-half_pi);

        if(!vec_cct->empty()) vec_cct->clear();
        vec_cct->push_back(this->CC_LN(q_s, *t));
        vec_cct->push_back(this->CC_RP(&(vec_cct->back()).q_g, *u));
        vec_cct->push_back(this->CC_LP(&(vec_cct->back()).q_g, *u));
        vec_cct->push_back(this->CC_RN(&(vec_cct->back()).q_g, *v));

        *err=this->distance(&(vec_cct->back()).q_g, q_g);
        return true;
    }
    else
        return false;
}
//27
bool CCStateSpace::RPLNRNLP(CCStateSpace::State *q_s,CCStateSpace::State *q_g, vector<CCStateSpace::CCTurn> *vec_cct, double *t, double *u, double *v, double *err, CCStateSpace::CCPathType *type)
{
    double l_max = 2*(1+2*cos(this->miu))*this->r;
    double l= dist2d(q_g->y_o_ln - q_s->y_o_rp, q_g->x_o_ln - q_s->x_o_rp);
    if (l <= l_max + CC_EPS)
    {
        *type=RpLnRnLp;
        double gamma=mod2pi2(atan2(q_g->y_o_ln - q_s->y_o_rp, q_g->x_o_ln - q_s->x_o_rp));
        double l1=l/4+this->r*this->r*(4*cos(this->miu)*cos(this->miu)-1)/l;
        double l2=l/2-l1;
        double beta1=(l>l_max-CC_EPS)?0:acos(l1/(2*this->r*cos(this->miu)));
        double beta2=(l>l_max-CC_EPS)?0:acos(l2/(this->r));

        *t=abs(mod2pi3(gamma-beta1-half_pi-q_s->theta));
        *u=mod2pi2(pi-beta1-beta2-this->miu);
        *v=mod2pi2(q_g->theta-gamma+beta1+half_pi);

        if(!vec_cct->empty()) vec_cct->clear();
        vec_cct->push_back(this->CC_RP(q_s, *t));
        vec_cct->push_back(this->CC_LN(&(vec_cct->back()).q_g, *u));
        vec_cct->push_back(this->CC_RN(&(vec_cct->back()).q_g, *u));
        vec_cct->push_back(this->CC_LP(&(vec_cct->back()).q_g, *v));

        *err=this->distance(&(vec_cct->back()).q_g, q_g);
        return true;
    }
    else
        return false;
}
//28
bool CCStateSpace::RNLPRPLN(CCStateSpace::State *q_s,CCStateSpace::State *q_g, vector<CCStateSpace::CCTurn> *vec_cct, double *t, double *u, double *v, double *err, CCStateSpace::CCPathType *type)
{
    double l_max = 2*(1+2*cos(this->miu))*this->r;
    double l= dist2d(q_g->y_o_lp - q_s->y_o_rn, q_g->x_o_lp - q_s->x_o_rn);
    if (l <= l_max + CC_EPS)
    {
        *type=RnLpRpLn;
        double gamma=mod2pi2(atan2(q_g->y_o_lp - q_s->y_o_rn, q_g->x_o_lp - q_s->x_o_rn));
        double l1=l/4+this->r*this->r*(4*cos(this->miu)*cos(this->miu)-1)/l;
        double l2=l/2-l1;
        double beta1=(l>l_max-CC_EPS)?0:acos(l1/(2*this->r*cos(this->miu)));
        double beta2=(l>l_max-CC_EPS)?0:acos(l2/(this->r));

        *t=mod2pi2(gamma+beta1-half_pi-q_s->theta);
        *u=mod2pi2(pi-beta1-beta2-this->miu);
        *v=abs(mod2pi3(q_g->theta-gamma-beta1+half_pi));

        if(!vec_cct->empty()) vec_cct->clear();
        vec_cct->push_back(this->CC_RN(q_s, *t));
        vec_cct->push_back(this->CC_LP(&(vec_cct->back()).q_g, *u));
        vec_cct->push_back(this->CC_RP(&(vec_cct->back()).q_g, *u));
        vec_cct->push_back(this->CC_LN(&(vec_cct->back()).q_g, *v));

        *err=this->distance(&(vec_cct->back()).q_g, q_g);
        return true;
    }
    else
        return false;
}
//CCSC,CSCC,16-------------------------------------------------------
//CCSC1
//29
bool CCStateSpace::LPRNSNLN(CCStateSpace::State *q_s,CCStateSpace::State *q_g, vector<CCStateSpace::CCTurn> *vec_cct, double *t, double *u, double *v, double *err, CCStateSpace::CCPathType *type)
{
    double uu = sqrt(dist2(q_g->y_o_lp-q_s->y_o_lp, q_g->x_o_lp-q_s->x_o_lp)-4*this->r*this->r*cos(this->miu)*cos(this->miu))-2*this->r*(sin(this->miu)+cos(this->miu));
    if (uu > -CC_EPS)
    {
        *u=uu;
        *type=LpRnSnLn;
        double gamma=mod2pi2(atan2(q_g->y_o_lp-q_s->y_o_lp, q_g->x_o_lp-q_s->x_o_lp));
        double beta = asin(2*this->r*cos(this->miu)/dist2d(q_g->y_o_lp-q_s->y_o_lp, q_g->x_o_lp-q_s->x_o_lp));
        *t=mod2pi2(gamma+beta+half_pi-q_s->theta);
        *v=abs(mod2pi3(q_g->theta-gamma-beta-pi));

        if(!vec_cct->empty()) vec_cct->clear();
        vec_cct->push_back(this->CC_LP(q_s, *t));
        vec_cct->push_back(this->CC_RN(&(vec_cct->back()).q_g, half_pi));
        if(*u>CC_EPS) vec_cct->push_back(this->CC_SN(&(vec_cct->back()).q_g, *u));
        vec_cct->push_back(this->CC_LN(&(vec_cct->back()).q_g, *v));

        *err=this->distance(&(vec_cct->back()).q_g, q_g);
        return true;
    }
    else
        return false;
}
//30
bool CCStateSpace::LNRPSPLP(CCStateSpace::State *q_s,CCStateSpace::State *q_g, vector<CCStateSpace::CCTurn> *vec_cct, double *t, double *u, double *v, double *err, CCStateSpace::CCPathType *type)
{
    double uu = sqrt(dist2(q_g->y_o_ln-q_s->y_o_ln, q_g->x_o_ln-q_s->x_o_ln)-4*this->r*this->r*cos(this->miu)*cos(this->miu))-2*this->r*(sin(this->miu)+cos(this->miu));
    if (uu > -CC_EPS)
    {
        *u=uu;
        *type=LnRpSpLp;
        double gamma=mod2pi2(atan2(q_g->y_o_ln-q_s->y_o_ln, q_g->x_o_ln-q_s->x_o_ln));
        double beta = asin(2*this->r*cos(this->miu)/dist2d(q_g->y_o_ln-q_s->y_o_ln, q_g->x_o_ln-q_s->x_o_ln));
        *t=abs(mod2pi3(gamma-beta+half_pi-q_s->theta));
        *v=mod2pi2(q_g->theta-gamma+beta);

        if(!vec_cct->empty()) vec_cct->clear();
        vec_cct->push_back(this->CC_LN(q_s, *t));
        vec_cct->push_back(this->CC_RP(&(vec_cct->back()).q_g, half_pi));
        if(*u>CC_EPS) vec_cct->push_back(this->CC_SP(&(vec_cct->back()).q_g, *u));
        vec_cct->push_back(this->CC_LP(&(vec_cct->back()).q_g, *v));

        *err=this->distance(&(vec_cct->back()).q_g, q_g);
        return true;
    }
    else
        return false;
}
//31
bool CCStateSpace::RPLNSNRN(CCStateSpace::State *q_s,CCStateSpace::State *q_g, vector<CCStateSpace::CCTurn> *vec_cct, double *t, double *u, double *v, double *err, CCStateSpace::CCPathType *type)
{
    double uu = sqrt(dist2(q_g->y_o_rp-q_s->y_o_rp, q_g->x_o_rp-q_s->x_o_rp)-4*this->r*this->r*cos(this->miu)*cos(this->miu))-2*this->r*(sin(this->miu)+cos(this->miu));
    if (uu > -CC_EPS)
    {
        *u=uu;
        *type=RpLnSnRn;
        double gamma=mod2pi2(atan2(q_g->y_o_rp-q_s->y_o_rp, q_g->x_o_rp-q_s->x_o_rp));
        double beta = asin(2*this->r*cos(this->miu)/dist2d(q_g->y_o_rp-q_s->y_o_rp, q_g->x_o_rp-q_s->x_o_rp));
        *t=abs(mod2pi3(gamma-beta-half_pi-q_s->theta));
        *v=mod2pi2(q_g->theta-gamma+beta+pi);

        if(!vec_cct->empty()) vec_cct->clear();
        vec_cct->push_back(this->CC_RP(q_s, *t));
        vec_cct->push_back(this->CC_LN(&(vec_cct->back()).q_g, half_pi));
        if(*u>CC_EPS) vec_cct->push_back(this->CC_SN(&(vec_cct->back()).q_g, *u));
        vec_cct->push_back(this->CC_RN(&(vec_cct->back()).q_g, *v));

        *err=this->distance(&(vec_cct->back()).q_g, q_g);
        return true;
    }
    else
        return false;
}
//32
bool CCStateSpace::RNLPSPRP(CCStateSpace::State *q_s,CCStateSpace::State *q_g, vector<CCStateSpace::CCTurn> *vec_cct, double *t, double *u, double *v, double *err, CCStateSpace::CCPathType *type)
{
    double uu = sqrt(dist2(q_g->y_o_rn-q_s->y_o_rn, q_g->x_o_rn-q_s->x_o_rn)-4*this->r*this->r*cos(this->miu)*cos(this->miu))-2*this->r*(sin(this->miu)+cos(this->miu));
    if (uu > -CC_EPS)
    {
        *u=uu;
        *type=RnLpSpRp;
        double gamma=mod2pi2(atan2(q_g->y_o_rn-q_s->y_o_rn, q_g->x_o_rn-q_s->x_o_rn));
        double beta = asin(2*this->r*cos(this->miu)/dist2d(q_g->y_o_rn-q_s->y_o_rn, q_g->x_o_rn-q_s->x_o_rn));
        *t=mod2pi2(gamma+beta-half_pi-q_s->theta);
        *v=abs(mod2pi3(q_g->theta-gamma-beta));

        if(!vec_cct->empty()) vec_cct->clear();
        vec_cct->push_back(this->CC_RN(q_s, *t));
        vec_cct->push_back(this->CC_LP(&(vec_cct->back()).q_g, half_pi));
        if(*u>CC_EPS) vec_cct->push_back(this->CC_SP(&(vec_cct->back()).q_g, *u));
        vec_cct->push_back(this->CC_RP(&(vec_cct->back()).q_g, *v));

        *err=this->distance(&(vec_cct->back()).q_g, q_g);
        return true;
    }
    else
        return false;
}
//CSCC1
//33
bool CCStateSpace::LNSNRNLP(CCStateSpace::State *q_s,CCStateSpace::State *q_g, vector<CCStateSpace::CCTurn> *vec_cct, double *t, double *u, double *v, double *err, CCStateSpace::CCPathType *type)
{
    double uu = sqrt(dist2(q_g->y_o_ln-q_s->y_o_ln, q_g->x_o_ln-q_s->x_o_ln)-4*this->r*this->r*cos(this->miu)*cos(this->miu))-2*this->r*(sin(this->miu)+cos(this->miu));
    if (uu > -CC_EPS)
    {
        *u=uu;
        *type=LnSnRnLp;
        double gamma=mod2pi2(atan2(q_g->y_o_ln-q_s->y_o_ln, q_g->x_o_ln-q_s->x_o_ln));
        double beta = asin(2*this->r*cos(this->miu)/dist2d(q_g->y_o_ln-q_s->y_o_ln, q_g->x_o_ln-q_s->x_o_ln));
        *t=abs(mod2pi3(gamma-beta+pi-q_s->theta));
        *v=mod2pi2(q_g->theta-gamma+beta-half_pi);

        if(!vec_cct->empty()) vec_cct->clear();
        vec_cct->push_back(this->CC_LN(q_s, *t));
        if(*u>CC_EPS) vec_cct->push_back(this->CC_SN(&(vec_cct->back()).q_g, *u));
        vec_cct->push_back(this->CC_RN(&(vec_cct->back()).q_g, half_pi));       
        vec_cct->push_back(this->CC_LP(&(vec_cct->back()).q_g, *v));

        *err=this->distance(&(vec_cct->back()).q_g, q_g);
        return true;
    }
    else
        return false;
}
//34
bool CCStateSpace::LPSPRPLN(CCStateSpace::State *q_s,CCStateSpace::State *q_g, vector<CCStateSpace::CCTurn> *vec_cct, double *t, double *u, double *v, double *err, CCStateSpace::CCPathType *type)
{
    double uu = sqrt(dist2(q_g->y_o_lp-q_s->y_o_lp, q_g->x_o_lp-q_s->x_o_lp)-4*this->r*this->r*cos(this->miu)*cos(this->miu))-2*this->r*(sin(this->miu)+cos(this->miu));
    if (uu > -CC_EPS)
    {
        *u=uu;
        *type=LpSpRpLn;
        double gamma=mod2pi2(atan2(q_g->y_o_lp-q_s->y_o_lp, q_g->x_o_lp-q_s->x_o_lp));
        double beta = asin(2*this->r*cos(this->miu)/dist2d(q_g->y_o_lp-q_s->y_o_lp, q_g->x_o_lp-q_s->x_o_lp));
        *t=mod2pi2(gamma+beta-q_s->theta);
        *v=abs(mod2pi3(q_g->theta-gamma-beta-half_pi));

        if(!vec_cct->empty()) vec_cct->clear();
        vec_cct->push_back(this->CC_LP(q_s, *t));
        if(*u>CC_EPS) vec_cct->push_back(this->CC_SP(&(vec_cct->back()).q_g, *u));
        vec_cct->push_back(this->CC_RP(&(vec_cct->back()).q_g, half_pi));       
        vec_cct->push_back(this->CC_LN(&(vec_cct->back()).q_g, *v));

        *err=this->distance(&(vec_cct->back()).q_g, q_g);
        return true;
    }
    else
        return false;
}
//35
bool CCStateSpace::RNSNLNRP(CCStateSpace::State *q_s,CCStateSpace::State *q_g, vector<CCStateSpace::CCTurn> *vec_cct, double *t, double *u, double *v, double *err, CCStateSpace::CCPathType *type)
{
    double uu = sqrt(dist2(q_g->y_o_rn-q_s->y_o_rn, q_g->x_o_rn-q_s->x_o_rn)-4*this->r*this->r*cos(this->miu)*cos(this->miu))-2*this->r*(sin(this->miu)+cos(this->miu));
    if (uu > -CC_EPS)
    {
        *u=uu;
        *type=RnSnLnRp;
        double gamma=mod2pi2(atan2(q_g->y_o_rn-q_s->y_o_rn, q_g->x_o_rn-q_s->x_o_rn));
        double beta = asin(2*this->r*cos(this->miu)/dist2d(q_g->y_o_rn-q_s->y_o_rn, q_g->x_o_rn-q_s->x_o_rn));
        *t=mod2pi2(gamma+beta+pi-q_s->theta);
        *v=abs(mod2pi3(q_g->theta-gamma-beta-half_pi));

        if(!vec_cct->empty()) vec_cct->clear();
        vec_cct->push_back(this->CC_RN(q_s, *t));
        if(*u>CC_EPS) vec_cct->push_back(this->CC_SN(&(vec_cct->back()).q_g, *u));
        vec_cct->push_back(this->CC_LN(&(vec_cct->back()).q_g, half_pi));       
        vec_cct->push_back(this->CC_RP(&(vec_cct->back()).q_g, *v));

        *err=this->distance(&(vec_cct->back()).q_g, q_g);
        return true;
    }
    else
        return false;
}
//36
bool CCStateSpace::RPSPLPRN(CCStateSpace::State *q_s,CCStateSpace::State *q_g, vector<CCStateSpace::CCTurn> *vec_cct, double *t, double *u, double *v, double *err, CCStateSpace::CCPathType *type)
{
    double uu = sqrt(dist2(q_g->y_o_rp-q_s->y_o_rp, q_g->x_o_rp-q_s->x_o_rp)-4*this->r*this->r*cos(this->miu)*cos(this->miu))-2*this->r*(sin(this->miu)+cos(this->miu));
    if (uu > -CC_EPS)
    {
        *u=uu;
        *type=RpSpLpRn;
        double gamma=mod2pi2(atan2(q_g->y_o_rp-q_s->y_o_rp, q_g->x_o_rp-q_s->x_o_rp));
        double beta = asin(2*this->r*cos(this->miu)/dist2d(q_g->y_o_rp-q_s->y_o_rp, q_g->x_o_rp-q_s->x_o_rp));
        *t=mod2pi2(gamma-beta-q_s->theta);
        *v=abs(mod2pi3(q_g->theta-gamma+beta+half_pi));

        if(!vec_cct->empty()) vec_cct->clear();
        vec_cct->push_back(this->CC_RP(q_s, *t));
        if(*u>CC_EPS) vec_cct->push_back(this->CC_SP(&(vec_cct->back()).q_g, *u));
        vec_cct->push_back(this->CC_LP(&(vec_cct->back()).q_g, half_pi));       
        vec_cct->push_back(this->CC_RN(&(vec_cct->back()).q_g, *v));

        *err=this->distance(&(vec_cct->back()).q_g, q_g);
        return true;
    }
    else
        return false;
}
//CCSC2
//37
bool CCStateSpace::LPRNSNRN(CCStateSpace::State *q_s,CCStateSpace::State *q_g, vector<CCStateSpace::CCTurn> *vec_cct, double *t, double *u, double *v, double *err, CCStateSpace::CCPathType *type)
{
    double uu = dist2d(q_g->y_o_rp-q_s->y_o_lp, q_g->x_o_rp-q_s->x_o_lp)-2*this->r*(sin(this->miu)+cos(this->miu));
    if (uu > -CC_EPS)
    {
        *u=uu;
        *type=LpRnSnRn;
        double gamma=mod2pi2(atan2(q_g->y_o_rp-q_s->y_o_lp, q_g->x_o_rp-q_s->x_o_lp));
        *t=mod2pi2(gamma+half_pi-q_s->theta);
        *v=mod2pi2(q_g->theta-gamma-pi);

        if(!vec_cct->empty()) vec_cct->clear();
        vec_cct->push_back(this->CC_LP(q_s, *t));
        vec_cct->push_back(this->CC_RN(&(vec_cct->back()).q_g, half_pi));
        if(*u>CC_EPS) vec_cct->push_back(this->CC_SN(&(vec_cct->back()).q_g, *u));
        vec_cct->push_back(this->CC_RN(&(vec_cct->back()).q_g, *v));

        *err=this->distance(&(vec_cct->back()).q_g, q_g);
        return true;
    }
    else
        return false;
}
//38
bool CCStateSpace::LNRPSPRP(CCStateSpace::State *q_s,CCStateSpace::State *q_g, vector<CCStateSpace::CCTurn> *vec_cct, double *t, double *u, double *v, double *err, CCStateSpace::CCPathType *type)
{
    double uu = dist2d(q_g->y_o_rn-q_s->y_o_ln, q_g->x_o_rn-q_s->x_o_ln)-2*this->r*(sin(this->miu)+cos(this->miu));
    if (uu > -CC_EPS)
    {
        *u=uu;
        *type=LnRpSpRp;
        double gamma=mod2pi2(atan2(q_g->y_o_rn-q_s->y_o_ln, q_g->x_o_rn-q_s->x_o_ln));
        *t=abs(mod2pi3(gamma+half_pi-q_s->theta));
        *v=abs(mod2pi3(q_g->theta-gamma));

        if(!vec_cct->empty()) vec_cct->clear();
        vec_cct->push_back(this->CC_LN(q_s, *t));
        vec_cct->push_back(this->CC_RP(&(vec_cct->back()).q_g, half_pi));
        if(*u>CC_EPS) vec_cct->push_back(this->CC_SP(&(vec_cct->back()).q_g, *u));
        vec_cct->push_back(this->CC_RP(&(vec_cct->back()).q_g, *v));

        *err=this->distance(&(vec_cct->back()).q_g, q_g);
        return true;
    }
    else
        return false;
}
//39
bool CCStateSpace::RPLNSNLN(CCStateSpace::State *q_s,CCStateSpace::State *q_g, vector<CCStateSpace::CCTurn> *vec_cct, double *t, double *u, double *v, double *err, CCStateSpace::CCPathType *type)
{
    double uu = dist2d(q_g->y_o_lp-q_s->y_o_rp, q_g->x_o_lp-q_s->x_o_rp)-2*this->r*(sin(this->miu)+cos(this->miu));
    if (uu > -CC_EPS)
    {
        *u=uu;
        *type=RpLnSnLn;
        double gamma=mod2pi2(atan2(q_g->y_o_lp-q_s->y_o_rp, q_g->x_o_lp-q_s->x_o_rp));
        *t=abs(mod2pi3(gamma-half_pi-q_s->theta));
        *v=abs(mod2pi3(q_g->theta-gamma+pi));

        if(!vec_cct->empty()) vec_cct->clear();
        vec_cct->push_back(this->CC_RP(q_s, *t));
        vec_cct->push_back(this->CC_LN(&(vec_cct->back()).q_g, half_pi));
        if(*u>CC_EPS) vec_cct->push_back(this->CC_SN(&(vec_cct->back()).q_g, *u));
        vec_cct->push_back(this->CC_LN(&(vec_cct->back()).q_g, *v));

        *err=this->distance(&(vec_cct->back()).q_g, q_g);
        return true;
    }
    else
        return false;
}
//40
bool CCStateSpace::RNLPSPLP(CCStateSpace::State *q_s,CCStateSpace::State *q_g, vector<CCStateSpace::CCTurn> *vec_cct, double *t, double *u, double *v, double *err, CCStateSpace::CCPathType *type)
{
    double uu = dist2d(q_g->y_o_ln-q_s->y_o_rn, q_g->x_o_ln-q_s->x_o_rn)-2*this->r*(sin(this->miu)+cos(this->miu));
    if (uu > -CC_EPS)
    {
        *u=uu;
        *type=RnLpSpLp;
        double gamma=mod2pi2(atan2(q_g->y_o_ln-q_s->y_o_rn, q_g->x_o_ln-q_s->x_o_rn));
        *t=mod2pi2(gamma-half_pi-q_s->theta);
        *v=mod2pi2(q_g->theta-gamma);

        if(!vec_cct->empty()) vec_cct->clear();
        vec_cct->push_back(this->CC_RN(q_s, *t));
        vec_cct->push_back(this->CC_LP(&(vec_cct->back()).q_g, half_pi));
        if(*u>CC_EPS) vec_cct->push_back(this->CC_SP(&(vec_cct->back()).q_g, *u));
        vec_cct->push_back(this->CC_LP(&(vec_cct->back()).q_g, *v));

        *err=this->distance(&(vec_cct->back()).q_g, q_g);
        return true;
    }
    else
        return false;
}
//CSCC2
//41
bool CCStateSpace::RNSNRNLP(CCStateSpace::State *q_s,CCStateSpace::State *q_g, vector<CCStateSpace::CCTurn> *vec_cct, double *t, double *u, double *v, double *err, CCStateSpace::CCPathType *type)
{
    double uu = dist2d(q_g->y_o_ln-q_s->y_o_rn, q_g->x_o_ln-q_s->x_o_rn)-2*this->r*(sin(this->miu)+cos(this->miu));
    if (uu > -CC_EPS)
    {
        *u=uu;
        *type=RnSnRnLp;
        double gamma=mod2pi2(atan2(q_g->y_o_ln-q_s->y_o_rn, q_g->x_o_ln-q_s->x_o_rn));
        *t=mod2pi2(gamma-pi-q_s->theta);
        *v=mod2pi2(q_g->theta-gamma+half_pi);

        if(!vec_cct->empty()) vec_cct->clear();
        vec_cct->push_back(this->CC_RN(q_s, *t));
        if(*u>CC_EPS) vec_cct->push_back(this->CC_SN(&(vec_cct->back()).q_g, *u));
        vec_cct->push_back(this->CC_RN(&(vec_cct->back()).q_g, half_pi));
        vec_cct->push_back(this->CC_LP(&(vec_cct->back()).q_g, *v));

        *err=this->distance(&(vec_cct->back()).q_g, q_g);
        return true;
    }
    else
        return false;
}
//42
bool CCStateSpace::RPSPRPLN(CCStateSpace::State *q_s,CCStateSpace::State *q_g, vector<CCStateSpace::CCTurn> *vec_cct, double *t, double *u, double *v, double *err, CCStateSpace::CCPathType *type)
{
    double uu = dist2d(q_g->y_o_lp-q_s->y_o_rp, q_g->x_o_lp-q_s->x_o_rp)-2*this->r*(sin(this->miu)+cos(this->miu));
    if (uu > -CC_EPS)
    {
        *u=uu;
        *type=RpSpRpLn;
        double gamma=mod2pi2(atan2(q_g->y_o_lp-q_s->y_o_rp, q_g->x_o_lp-q_s->x_o_rp));
        *t=abs(mod2pi3(gamma-q_s->theta));
        *v=abs(mod2pi3(q_g->theta-gamma+half_pi));

        if(!vec_cct->empty()) vec_cct->clear();
        vec_cct->push_back(this->CC_RP(q_s, *t));
        if(*u>CC_EPS) vec_cct->push_back(this->CC_SP(&(vec_cct->back()).q_g, *u));
        vec_cct->push_back(this->CC_RP(&(vec_cct->back()).q_g, half_pi));
        vec_cct->push_back(this->CC_LN(&(vec_cct->back()).q_g, *v));

        *err=this->distance(&(vec_cct->back()).q_g, q_g);
        return true;
    }
    else
        return false;
}
//43
bool CCStateSpace::LNSNLNRP(CCStateSpace::State *q_s,CCStateSpace::State *q_g, vector<CCStateSpace::CCTurn> *vec_cct, double *t, double *u, double *v, double *err, CCStateSpace::CCPathType *type)
{
    double uu = dist2d(q_g->y_o_rn-q_s->y_o_ln, q_g->x_o_rn-q_s->x_o_ln)-2*this->r*(sin(this->miu)+cos(this->miu));
    if (uu > -CC_EPS)
    {
        *u=uu;
        *type=LnSnLnRp;
        double gamma=mod2pi2(atan2(q_g->y_o_rn-q_s->y_o_ln, q_g->x_o_rn-q_s->x_o_ln));
        *t=abs(mod2pi3(gamma+pi-q_s->theta));
        *v=abs(mod2pi3(q_g->theta-gamma-half_pi));

        if(!vec_cct->empty()) vec_cct->clear();
        vec_cct->push_back(this->CC_LN(q_s, *t));
        if(*u>CC_EPS) vec_cct->push_back(this->CC_SN(&(vec_cct->back()).q_g, *u));
        vec_cct->push_back(this->CC_LN(&(vec_cct->back()).q_g, half_pi));
        vec_cct->push_back(this->CC_RP(&(vec_cct->back()).q_g, *v));

        *err=this->distance(&(vec_cct->back()).q_g, q_g);
        return true;
    }
    else
        return false;
}
//44
bool CCStateSpace::LPSPLPRN(CCStateSpace::State *q_s,CCStateSpace::State *q_g, vector<CCStateSpace::CCTurn> *vec_cct, double *t, double *u, double *v, double *err, CCStateSpace::CCPathType *type)
{
    double uu = dist2d(q_g->y_o_rp-q_s->y_o_lp, q_g->x_o_rp-q_s->x_o_lp)-2*this->r*(sin(this->miu)+cos(this->miu));
    if (uu > -CC_EPS)
    {
        *u=uu;
        *type=LpSpLpRn;
        double gamma=mod2pi2(atan2(q_g->y_o_rp-q_s->y_o_lp, q_g->x_o_rp-q_s->x_o_lp));
        *t=mod2pi2(gamma-q_s->theta);
        *v=mod2pi2(q_g->theta-gamma-half_pi);

        if(!vec_cct->empty()) vec_cct->clear();
        vec_cct->push_back(this->CC_LP(q_s, *t));
        if(*u>CC_EPS) vec_cct->push_back(this->CC_SP(&(vec_cct->back()).q_g, *u));
        vec_cct->push_back(this->CC_LP(&(vec_cct->back()).q_g, half_pi));
        vec_cct->push_back(this->CC_RN(&(vec_cct->back()).q_g, *v));

        *err=this->distance(&(vec_cct->back()).q_g, q_g);
        return true;
    }
    else
        return false;
}
//CCSCC,4-------------------------------------------------------
//45
bool CCStateSpace::LPRNSNLNRP(CCStateSpace::State *q_s,CCStateSpace::State *q_g, vector<CCStateSpace::CCTurn> *vec_cct, double *t, double *u, double *v, double *err, CCStateSpace::CCPathType *type)
{
    double uu = sqrt(dist2(q_g->y_o_rn-q_s->y_o_lp, q_g->x_o_rn-q_s->x_o_lp)-4*this->r*this->r*cos(this->miu)*cos(this->miu))-2*this->r*(2*cos(this->miu)+sin(this->miu));
    if (uu > -CC_EPS)
    {
        *u=uu;
        *type=LpRnSnLnRp;
        double gamma=mod2pi2(atan2(q_g->y_o_rn-q_s->y_o_lp, q_g->x_o_rn-q_s->x_o_lp));
        double beta=mod2pi2(asin(2*this->r*cos(this->miu)/dist2d(q_g->y_o_rn-q_s->y_o_lp, q_g->x_o_rn-q_s->x_o_lp)));
        *t=mod2pi2(gamma+beta+half_pi-q_s->theta);
        *v=abs(mod2pi3(q_g->theta-gamma-beta-half_pi));

        if(!vec_cct->empty()) vec_cct->clear();
        vec_cct->push_back(this->CC_LP(q_s, *t));
        vec_cct->push_back(this->CC_RN(&(vec_cct->back()).q_g, half_pi));
        if(*u>CC_EPS) vec_cct->push_back(this->CC_SN(&(vec_cct->back()).q_g, *u));
        vec_cct->push_back(this->CC_LN(&(vec_cct->back()).q_g, half_pi));
        vec_cct->push_back(this->CC_RP(&(vec_cct->back()).q_g, *v));

        *err=this->distance(&(vec_cct->back()).q_g, q_g);
        return true;
    }
    else
        return false;
}
//46
bool CCStateSpace::LNRPSPLPRN(CCStateSpace::State *q_s,CCStateSpace::State *q_g, vector<CCStateSpace::CCTurn> *vec_cct, double *t, double *u, double *v, double *err, CCStateSpace::CCPathType *type)
{
    double uu = sqrt(dist2(q_g->y_o_rp-q_s->y_o_ln, q_g->x_o_rp-q_s->x_o_ln)-4*this->r*this->r*cos(this->miu)*cos(this->miu))-2*this->r*(2*cos(this->miu)+sin(this->miu));
    if (uu > -CC_EPS)
    {
        *u=uu;
        *type=LnRpSpLpRn;
        double gamma=mod2pi2(atan2(q_g->y_o_rp-q_s->y_o_ln, q_g->x_o_rp-q_s->x_o_ln));
        double beta=mod2pi2(asin(2*this->r*cos(this->miu)/dist2d(q_g->y_o_rp-q_s->y_o_ln, q_g->x_o_rp-q_s->x_o_ln)));
        *t=abs(mod2pi3(gamma-beta+half_pi-q_s->theta));
        *v=mod2pi2(q_g->theta-gamma+beta-half_pi);

        if(!vec_cct->empty()) vec_cct->clear();
        vec_cct->push_back(this->CC_LN(q_s, *t));
        vec_cct->push_back(this->CC_RP(&(vec_cct->back()).q_g, half_pi));
        if(*u>CC_EPS) vec_cct->push_back(this->CC_SP(&(vec_cct->back()).q_g, *u));
        vec_cct->push_back(this->CC_LP(&(vec_cct->back()).q_g, half_pi));
        vec_cct->push_back(this->CC_RN(&(vec_cct->back()).q_g, *v));

        *err=this->distance(&(vec_cct->back()).q_g, q_g);
        return true;
    }
    else
        return false;
}
//47
bool CCStateSpace::RPLNSNRNLP(CCStateSpace::State *q_s,CCStateSpace::State *q_g, vector<CCStateSpace::CCTurn> *vec_cct, double *t, double *u, double *v, double *err, CCStateSpace::CCPathType *type)
{
    double uu = sqrt(dist2(q_g->y_o_ln-q_s->y_o_rp, q_g->x_o_ln-q_s->x_o_rp)-4*this->r*this->r*cos(this->miu)*cos(this->miu))-2*this->r*(2*cos(this->miu)+sin(this->miu));
    if (uu > -CC_EPS)
    {
        *u=uu;
        *type=RpLnSnRnLp;
        double gamma=mod2pi2(atan2(q_g->y_o_ln-q_s->y_o_rp, q_g->x_o_ln-q_s->x_o_rp));
        double beta=mod2pi2(asin(2*this->r*cos(this->miu)/dist2d(q_g->y_o_ln-q_s->y_o_rp, q_g->x_o_ln-q_s->x_o_rp)));
        *t=abs(mod2pi3(gamma-beta-half_pi-q_s->theta));
        *v=mod2pi2(q_g->theta-gamma+beta+half_pi);

        if(!vec_cct->empty()) vec_cct->clear();
        vec_cct->push_back(this->CC_RP(q_s, *t));
        vec_cct->push_back(this->CC_LN(&(vec_cct->back()).q_g, half_pi));
        if(*u>CC_EPS) vec_cct->push_back(this->CC_SN(&(vec_cct->back()).q_g, *u));
        vec_cct->push_back(this->CC_RN(&(vec_cct->back()).q_g, half_pi));
        vec_cct->push_back(this->CC_LP(&(vec_cct->back()).q_g, *v));

        *err=this->distance(&(vec_cct->back()).q_g, q_g);
        return true;
    }
    else
        return false;
}
//48
bool CCStateSpace::RNLPSPRPLN(CCStateSpace::State *q_s,CCStateSpace::State *q_g, vector<CCStateSpace::CCTurn> *vec_cct, double *t, double *u, double *v, double *err, CCStateSpace::CCPathType *type)
{
    double uu = sqrt(dist2(q_g->y_o_lp-q_s->y_o_rn, q_g->x_o_lp-q_s->x_o_rn)-4*this->r*this->r*cos(this->miu)*cos(this->miu))-2*this->r*(2*cos(this->miu)+sin(this->miu));
    if (uu > -CC_EPS)
    {
        *u=uu;
        *type=RnLpSpRpLn;
        double gamma=mod2pi2(atan2(q_g->y_o_lp-q_s->y_o_rn, q_g->x_o_lp-q_s->x_o_rn));
        double beta=mod2pi2(asin(2*this->r*cos(this->miu)/dist2d(q_g->y_o_lp-q_s->y_o_rn, q_g->x_o_lp-q_s->x_o_rn)));
        *t=mod2pi2(gamma+beta-half_pi-q_s->theta);
        *v=abs(mod2pi3(q_g->theta-gamma-beta+half_pi));

        if(!vec_cct->empty()) vec_cct->clear();
        vec_cct->push_back(this->CC_RN(q_s, *t));
        vec_cct->push_back(this->CC_LP(&(vec_cct->back()).q_g, half_pi));
        if(*u>CC_EPS) vec_cct->push_back(this->CC_SP(&(vec_cct->back()).q_g, *u));
        vec_cct->push_back(this->CC_RP(&(vec_cct->back()).q_g, half_pi));
        vec_cct->push_back(this->CC_LN(&(vec_cct->back()).q_g, *v));

        *err=this->distance(&(vec_cct->back()).q_g, q_g);
        return true;
    }
    else
        return false;
}
//----------------------------------------------------------------
inline double len_cct(vector<CCStateSpace::CCTurn> *vec)
{
    double l = 0;
    if (!vec->empty())
    {
        for (size_t i = 0; i<vec->size(); ++i)
            l += (*vec)[i].length;
    }
    return l;
}
//Topology-------------------------------------------------------
//distance between q_s and q_g little than half_pi/space.k_max, to keep q_s_ref.theta in (-pi/2, pi/2)
void CCStateSpace::TOPOLOGY(CCStateSpace::State *q_s,CCStateSpace::State *q_g, vector<vector<CCStateSpace::CCTurn>> *vec, vector<double> *vt, vector<double> *vu, vector<double> *vv, vector<double> *vlength, vector<double> *verr, vector<CCStateSpace::CCPathType> *vtype)
{
    if(this->distance(q_s,q_g)<two_pi/this->k_max &&abs(mod2pi(q_g->theta-q_s->theta))<half_pi)
    {
        vector<CCStateSpace::CCTurn> vec_cct;
        double t,u,v;

        double d_m=2*root_pi/sqrt(this->sigma_max)*D(this->theta_i)*sin(this->theta_i)/cos(this->delta_min);
        CCStateSpace::State q_s_ref(q_g->local_to_ref(q_s, *this));
        double abs_theta_s=abs(q_s_ref.theta);
        double d_tmp=D(abs_theta_s/2)/this->r;
        if(abs_theta_s>this->delta_min-CC_EPS)
            t=4*pi*d_tmp*d_tmp;
        else
            t=this->sigma_max;

        //1,two clothoids
        if(q_s_ref.theta > CC_EPS)
        {
            vec_cct.push_back(CCStateSpace::CCTurn(CCStateSpace::ccTurnType[10], *q_s, abs_theta_s, t, *this)); //LN
        }
        else if(q_s_ref.theta < -CC_EPS)
        {
            vec_cct.push_back(CCStateSpace::CCTurn(CCStateSpace::ccTurnType[13], *q_s, abs_theta_s, t, *this)); //RN
        }
        //2, line
        CCStateSpace::State *q=q_s;
        CCStateSpace::State *q_r=&q_s_ref;
        if(!vec_cct.empty()) 
            {
                q=&(vec_cct.back().q_g);
                *q_r=q_g->local_to_ref(q, *this);
            }
        if(q_r->x < -CC_EPS)
        {
            vec_cct.push_back(CC_SP(q,-q_r->x));
        }
        else if(q_r->x > CC_EPS)
        {
            vec_cct.push_back(CC_SN(q, q_r->x));
        }
        //3, two clothoids + line + two clothoids
        if(!vec_cct.empty()) 
            {
                q=&(vec_cct.back().q_g);
                *q_r=q_g->local_to_ref(q, *this);
            }
        u=solve_alpha(abs(q_r->y)/2, d_m, *this);
        double l=abs(q_r->y)/tan(u);
        if(u > this->theta_i)
        {
            d_tmp=D((u))/this->r;
            v=4*pi*d_tmp*d_tmp;
        }
        else
        {
            v=this->sigma_max;
        }
        if(q_r->y > CC_EPS)
        {
            //Lp
            vec_cct.push_back(CCStateSpace::CCTurn(CCStateSpace::ccTurnType[4], *q, 2*(u), v, *this));
            //Sn
            vec_cct.push_back(CCStateSpace::CCTurn(CCStateSpace::ccTurnType[1], (vec_cct.back()).q_g, l,*this));
            //Rp
            vec_cct.push_back(CCStateSpace::CCTurn(CCStateSpace::ccTurnType[7], (vec_cct.back()).q_g, 2*(u), v, *this));
        }
        else if(q_r->y < -CC_EPS)
        {
            //Rp
            vec_cct.push_back(CCStateSpace::CCTurn(CCStateSpace::ccTurnType[7], *q, 2*(u), v, *this));
            //Sn
            vec_cct.push_back(CCStateSpace::CCTurn(CCStateSpace::ccTurnType[1], (vec_cct.back()).q_g, l,*this));
            //Lp
            vec_cct.push_back(CCStateSpace::CCTurn(CCStateSpace::ccTurnType[4], (vec_cct.back()).q_g, 2*(u), v, *this));
        }

        if(!vec_cct.empty())
        {
            vec->push_back(vec_cct);
            vt->push_back(t); vu->push_back(u); vv->push_back(v);
            vlength->push_back(len_cct(&vec_cct));
            verr->push_back(this->distance(&(vec_cct.back()).q_g, q_g));
            vtype->push_back(Topology);
        }           
    }
}
//----------------------CSC\CCC\CCCC\CCSC\CCSCC-----------------------------------------------

void CCStateSpace::CSC(CCStateSpace::State *q_s,CCStateSpace::State *q_g, vector<vector<CCStateSpace::CCTurn>> *vec_cct, vector<double> *t, vector<double> *u, vector<double> *v, vector<double> *length, vector<double> *err, vector<CCStateSpace::CCPathType> *type)
{
    vector<CCStateSpace::CCTurn> vec_tmp;
    double t1,u1,v1,err1;
    CCStateSpace::CCPathType tp;
    if(LPSPLP(q_s,q_g,&vec_tmp,&t1,&u1,&v1,&err1,&tp))
    {
        vec_cct->push_back(vec_tmp);
        t->push_back(t1); u->push_back(u1); v->push_back(v1);
        length->push_back(len_cct(&vec_tmp));
        err->push_back(err1);
        type->push_back(tp);
    }
    else if(RPSPRP(q_s,q_g,&vec_tmp,&t1,&u1,&v1,&err1,&tp))
    {
        vec_cct->push_back(vec_tmp);
        t->push_back(t1); u->push_back(u1); v->push_back(v1);
        length->push_back(len_cct(&vec_tmp));
        err->push_back(err1);
        type->push_back(tp);
    }
    else if(LNSNLN(q_s,q_g,&vec_tmp,&t1,&u1,&v1,&err1,&tp))
    {
        vec_cct->push_back(vec_tmp);
        t->push_back(t1); u->push_back(u1); v->push_back(v1);
        length->push_back(len_cct(&vec_tmp));
        err->push_back(err1);
        type->push_back(tp);
    }
    else if(RNSNRN(q_s,q_g,&vec_tmp,&t1,&u1,&v1,&err1,&tp))
    {
        vec_cct->push_back(vec_tmp);
        t->push_back(t1); u->push_back(u1); v->push_back(v1);
        length->push_back(len_cct(&vec_tmp));
        err->push_back(err1);
        type->push_back(tp);
    }
    else if(LPSPRP(q_s,q_g,&vec_tmp,&t1,&u1,&v1,&err1,&tp))
    {
        vec_cct->push_back(vec_tmp);
        t->push_back(t1); u->push_back(u1); v->push_back(v1);
        length->push_back(len_cct(&vec_tmp));
        err->push_back(err1);
        type->push_back(tp);
    }
    else if(RPSPLP(q_s,q_g,&vec_tmp,&t1,&u1,&v1,&err1,&tp))
    {
        vec_cct->push_back(vec_tmp);
        t->push_back(t1); u->push_back(u1); v->push_back(v1);
        length->push_back(len_cct(&vec_tmp));
        err->push_back(err1);
        type->push_back(tp);
    }
    else if(LNSNRN(q_s,q_g,&vec_tmp,&t1,&u1,&v1,&err1,&tp))
    {
        vec_cct->push_back(vec_tmp);
        t->push_back(t1); u->push_back(u1); v->push_back(v1);
        length->push_back(len_cct(&vec_tmp));
        err->push_back(err1);
        type->push_back(tp);
    }
    else if(RNSNLN(q_s,q_g,&vec_tmp,&t1,&u1,&v1,&err1,&tp))
    {
        vec_cct->push_back(vec_tmp);
        t->push_back(t1); u->push_back(u1); v->push_back(v1);
        length->push_back(len_cct(&vec_tmp));
        err->push_back(err1);
        type->push_back(tp);
    }
}
void CCStateSpace::CCC(CCStateSpace::State *q_s, CCStateSpace::State *q_g, vector<vector<CCStateSpace::CCTurn>> *vec_cct, vector<double> *t, vector<double> *u, vector<double> *v, vector<double> *length, vector<double> *err, vector<CCStateSpace::CCPathType> *type)
{
    vector<CCStateSpace::CCTurn> vec_tmp;
    double t1, u1, v1, err1;
    CCStateSpace::CCPathType tp;
    if (LPRNLP(q_s, q_g, &vec_tmp, &t1, &u1, &v1, &err1, &tp))
    {
        vec_cct->push_back(vec_tmp);
        t->push_back(t1); u->push_back(u1); v->push_back(v1);
        length->push_back(len_cct(&vec_tmp));
        err->push_back(err1);
        type->push_back(tp);
    }
    else if (RPLNRP(q_s, q_g, &vec_tmp, &t1, &u1, &v1, &err1, &tp))
    {
        vec_cct->push_back(vec_tmp);
        t->push_back(t1); u->push_back(u1); v->push_back(v1);
        length->push_back(len_cct(&vec_tmp));
        err->push_back(err1);
        type->push_back(tp);
    }
    else if (LNRPLN(q_s, q_g, &vec_tmp, &t1, &u1, &v1, &err1, &tp))
    {
        vec_cct->push_back(vec_tmp);
        t->push_back(t1); u->push_back(u1); v->push_back(v1);
        length->push_back(len_cct(&vec_tmp));
        err->push_back(err1);
        type->push_back(tp);
    }
    else if (RNLPRN(q_s, q_g, &vec_tmp, &t1, &u1, &v1, &err1, &tp))
    {
        vec_cct->push_back(vec_tmp);
        t->push_back(t1); u->push_back(u1); v->push_back(v1);
        length->push_back(len_cct(&vec_tmp));
        err->push_back(err1);
        type->push_back(tp);
    }
    else if (LPRNLN(q_s, q_g, &vec_tmp, &t1, &u1, &v1, &err1, &tp))
    {
        vec_cct->push_back(vec_tmp);
        t->push_back(t1); u->push_back(u1); v->push_back(v1);
        length->push_back(len_cct(&vec_tmp));
        err->push_back(err1);
        type->push_back(tp);
    }
    else if (LNRPLP(q_s, q_g, &vec_tmp, &t1, &u1, &v1, &err1, &tp))
    {
        vec_cct->push_back(vec_tmp);
        t->push_back(t1); u->push_back(u1); v->push_back(v1);
        length->push_back(len_cct(&vec_tmp));
        err->push_back(err1);
        type->push_back(tp);
    }
    else if (RPLNRN(q_s, q_g, &vec_tmp, &t1, &u1, &v1, &err1, &tp))
    {
        vec_cct->push_back(vec_tmp);
        t->push_back(t1); u->push_back(u1); v->push_back(v1);
        length->push_back(len_cct(&vec_tmp));
        err->push_back(err1);
        type->push_back(tp);
    }
    else if (RNLPRP(q_s, q_g, &vec_tmp, &t1, &u1, &v1, &err1, &tp))
    {
        vec_cct->push_back(vec_tmp);
        t->push_back(t1); u->push_back(u1); v->push_back(v1);
        length->push_back(len_cct(&vec_tmp));
        err->push_back(err1);
        type->push_back(tp);
    }
    else if (LNRNLP(q_s, q_g, &vec_tmp, &t1, &u1, &v1, &err1, &tp))
    {
        vec_cct->push_back(vec_tmp);
        t->push_back(t1); u->push_back(u1); v->push_back(v1);
        length->push_back(len_cct(&vec_tmp));
        err->push_back(err1);
        type->push_back(tp);
    }
    else if (LPRPLN(q_s, q_g, &vec_tmp, &t1, &u1, &v1, &err1, &tp))
    {
        vec_cct->push_back(vec_tmp);
        t->push_back(t1); u->push_back(u1); v->push_back(v1);
        length->push_back(len_cct(&vec_tmp));
        err->push_back(err1);
        type->push_back(tp);
    }
    else if (RNLNRP(q_s, q_g, &vec_tmp, &t1, &u1, &v1, &err1, &tp))
    {
        vec_cct->push_back(vec_tmp);
        t->push_back(t1); u->push_back(u1); v->push_back(v1);
        length->push_back(len_cct(&vec_tmp));
        err->push_back(err1);
        type->push_back(tp);
    }
    else if (RPLPRN(q_s, q_g, &vec_tmp, &t1, &u1, &v1, &err1, &tp))
    {
        vec_cct->push_back(vec_tmp);
        t->push_back(t1); u->push_back(u1); v->push_back(v1);
        length->push_back(len_cct(&vec_tmp));
        err->push_back(err1);
        type->push_back(tp);
    }

}
void CCStateSpace::CCCC(CCStateSpace::State *q_s,CCStateSpace::State *q_g, vector<vector<CCStateSpace::CCTurn>> *vec_cct, vector<double> *t, vector<double> *u, vector<double> *v, vector<double> *length, vector<double> *err, vector<CCStateSpace::CCPathType> *type)
{
    vector<CCStateSpace::CCTurn> vec_tmp;
    double t1,u1,v1,err1;
    CCStateSpace::CCPathType tp;
    if(LPRPLNRN(q_s,q_g,&vec_tmp,&t1,&u1,&v1,&err1,&tp))
    {
        vec_cct->push_back(vec_tmp);
        t->push_back(t1); u->push_back(u1); v->push_back(v1);
        length->push_back(len_cct(&vec_tmp));
        err->push_back(err1);
        type->push_back(tp);
    }
    else if(LNRNLPRP(q_s,q_g,&vec_tmp,&t1,&u1,&v1,&err1,&tp))
    {
        vec_cct->push_back(vec_tmp);
        t->push_back(t1); u->push_back(u1); v->push_back(v1);
        length->push_back(len_cct(&vec_tmp));
        err->push_back(err1);
        type->push_back(tp);
    }
    else if(RPLPRNLN(q_s,q_g,&vec_tmp,&t1,&u1,&v1,&err1,&tp))
    {
        vec_cct->push_back(vec_tmp);
        t->push_back(t1); u->push_back(u1); v->push_back(v1);
        length->push_back(len_cct(&vec_tmp));
        err->push_back(err1);
        type->push_back(tp);
    }
    else if(RNLNRPLP(q_s,q_g,&vec_tmp,&t1,&u1,&v1,&err1,&tp))
    {
        vec_cct->push_back(vec_tmp);
        t->push_back(t1); u->push_back(u1); v->push_back(v1);
        length->push_back(len_cct(&vec_tmp));
        err->push_back(err1);
        type->push_back(tp);
    }
    else if(LPRNLNRP(q_s,q_g,&vec_tmp,&t1,&u1,&v1,&err1,&tp))
    {
        vec_cct->push_back(vec_tmp);
        t->push_back(t1); u->push_back(u1); v->push_back(v1);
        length->push_back(len_cct(&vec_tmp));
        err->push_back(err1);
        type->push_back(tp);
    }
    else if(LNRPLPRN(q_s,q_g,&vec_tmp,&t1,&u1,&v1,&err1,&tp))
    {
        vec_cct->push_back(vec_tmp);
        t->push_back(t1); u->push_back(u1); v->push_back(v1);
        length->push_back(len_cct(&vec_tmp));
        err->push_back(err1);
        type->push_back(tp);
    }
    else if(RPLNRNLP(q_s,q_g,&vec_tmp,&t1,&u1,&v1,&err1,&tp))
    {
        vec_cct->push_back(vec_tmp);
        t->push_back(t1); u->push_back(u1); v->push_back(v1);
        length->push_back(len_cct(&vec_tmp));
        err->push_back(err1);
        type->push_back(tp);
    }
    else if(RNLPRPLN(q_s,q_g,&vec_tmp,&t1,&u1,&v1,&err1,&tp))
    {
        vec_cct->push_back(vec_tmp);
        t->push_back(t1); u->push_back(u1); v->push_back(v1);
        length->push_back(len_cct(&vec_tmp));
        err->push_back(err1);
        type->push_back(tp);
    }
}
void CCStateSpace::CCSC(CCStateSpace::State *q_s,CCStateSpace::State *q_g, vector<vector<CCStateSpace::CCTurn>> *vec_cct, vector<double> *t, vector<double> *u, vector<double> *v, vector<double> *length, vector<double> *err, vector<CCStateSpace::CCPathType> *type)
{
    vector<CCStateSpace::CCTurn> vec_tmp;
    double t1,u1,v1,err1;
    CCStateSpace::CCPathType tp;
    if(LPRNSNLN(q_s,q_g,&vec_tmp,&t1,&u1,&v1,&err1,&tp))
    {
        vec_cct->push_back(vec_tmp);
        t->push_back(t1); u->push_back(u1); v->push_back(v1);
        length->push_back(len_cct(&vec_tmp));
        err->push_back(err1);
        type->push_back(tp);
    }
    else if(LNRPSPLP(q_s,q_g,&vec_tmp,&t1,&u1,&v1,&err1,&tp))
    {
        vec_cct->push_back(vec_tmp);
        t->push_back(t1); u->push_back(u1); v->push_back(v1);
        length->push_back(len_cct(&vec_tmp));
        err->push_back(err1);
        type->push_back(tp);
    }
    else if(RPLNSNRN(q_s,q_g,&vec_tmp,&t1,&u1,&v1,&err1,&tp))
    {
        vec_cct->push_back(vec_tmp);
        t->push_back(t1); u->push_back(u1); v->push_back(v1);
        length->push_back(len_cct(&vec_tmp));
        err->push_back(err1);
        type->push_back(tp);
    }
    else if(RNLPSPRP(q_s,q_g,&vec_tmp,&t1,&u1,&v1,&err1,&tp))
    {
        vec_cct->push_back(vec_tmp);
        t->push_back(t1); u->push_back(u1); v->push_back(v1);
        length->push_back(len_cct(&vec_tmp));
        err->push_back(err1);
        type->push_back(tp);
    }
    else if(LNSNRNLP(q_s,q_g,&vec_tmp,&t1,&u1,&v1,&err1,&tp))
    {
        vec_cct->push_back(vec_tmp);
        t->push_back(t1); u->push_back(u1); v->push_back(v1);
        length->push_back(len_cct(&vec_tmp));
        err->push_back(err1);
        type->push_back(tp);
    }
    else if(LPSPRPLN(q_s,q_g,&vec_tmp,&t1,&u1,&v1,&err1,&tp))
    {
        vec_cct->push_back(vec_tmp);
        t->push_back(t1); u->push_back(u1); v->push_back(v1);
        length->push_back(len_cct(&vec_tmp));
        err->push_back(err1);
        type->push_back(tp);
    }
    else if(RNSNLNRP(q_s,q_g,&vec_tmp,&t1,&u1,&v1,&err1,&tp))
    {
        vec_cct->push_back(vec_tmp);
        t->push_back(t1); u->push_back(u1); v->push_back(v1);
        length->push_back(len_cct(&vec_tmp));
        err->push_back(err1);
        type->push_back(tp);
    }
    else if(RPSPLPRN(q_s,q_g,&vec_tmp,&t1,&u1,&v1,&err1,&tp))
    {
        vec_cct->push_back(vec_tmp);
        t->push_back(t1); u->push_back(u1); v->push_back(v1);
        length->push_back(len_cct(&vec_tmp));
        err->push_back(err1);
        type->push_back(tp);
    }
    else if(LPRNSNRN(q_s,q_g,&vec_tmp,&t1,&u1,&v1,&err1,&tp))
    {
        vec_cct->push_back(vec_tmp);
        t->push_back(t1); u->push_back(u1); v->push_back(v1);
        length->push_back(len_cct(&vec_tmp));
        err->push_back(err1);
        type->push_back(tp);
    }
    else if(LNRPSPRP(q_s,q_g,&vec_tmp,&t1,&u1,&v1,&err1,&tp))
    {
        vec_cct->push_back(vec_tmp);
        t->push_back(t1); u->push_back(u1); v->push_back(v1);
        length->push_back(len_cct(&vec_tmp));
        err->push_back(err1);
        type->push_back(tp);
    }
    else if(RPLNSNLN(q_s,q_g,&vec_tmp,&t1,&u1,&v1,&err1,&tp))
    {
        vec_cct->push_back(vec_tmp);
        t->push_back(t1); u->push_back(u1); v->push_back(v1);
        length->push_back(len_cct(&vec_tmp));
        err->push_back(err1);
        type->push_back(tp);
    }
    else if(RNLPSPLP(q_s,q_g,&vec_tmp,&t1,&u1,&v1,&err1,&tp))
    {
        vec_cct->push_back(vec_tmp);
        t->push_back(t1); u->push_back(u1); v->push_back(v1);
        length->push_back(len_cct(&vec_tmp));
        err->push_back(err1);
        type->push_back(tp);
    }
    else if(RNSNRNLP(q_s,q_g,&vec_tmp,&t1,&u1,&v1,&err1,&tp))
    {
        vec_cct->push_back(vec_tmp);
        t->push_back(t1); u->push_back(u1); v->push_back(v1);
        length->push_back(len_cct(&vec_tmp));
        err->push_back(err1);
        type->push_back(tp);
    }
    else if(RPSPRPLN(q_s,q_g,&vec_tmp,&t1,&u1,&v1,&err1,&tp))
    {
        vec_cct->push_back(vec_tmp);
        t->push_back(t1); u->push_back(u1); v->push_back(v1);
        length->push_back(len_cct(&vec_tmp));
        err->push_back(err1);
        type->push_back(tp);
    }
    else if(LNSNLNRP(q_s,q_g,&vec_tmp,&t1,&u1,&v1,&err1,&tp))
    {
        vec_cct->push_back(vec_tmp);
        t->push_back(t1); u->push_back(u1); v->push_back(v1);
        length->push_back(len_cct(&vec_tmp));
        err->push_back(err1);
        type->push_back(tp);
    }
    else if(LPSPLPRN(q_s,q_g,&vec_tmp,&t1,&u1,&v1,&err1,&tp))
    {
        vec_cct->push_back(vec_tmp);
        t->push_back(t1); u->push_back(u1); v->push_back(v1);
        length->push_back(len_cct(&vec_tmp));
        err->push_back(err1);
        type->push_back(tp);
    }
}

void CCStateSpace::CCSCC(CCStateSpace::State *q_s,CCStateSpace::State *q_g, vector<vector<CCStateSpace::CCTurn>> *vec_cct, vector<double> *t, vector<double> *u, vector<double> *v, vector<double> *length, vector<double> *err, vector<CCStateSpace::CCPathType> *type)
{
    vector<CCStateSpace::CCTurn> vec_tmp;
    double t1,u1,v1,err1;
    CCStateSpace::CCPathType tp;
    if(LPRNSNLNRP(q_s,q_g,&vec_tmp,&t1,&u1,&v1,&err1,&tp))
    {
        vec_cct->push_back(vec_tmp);
        t->push_back(t1); u->push_back(u1); v->push_back(v1);
        length->push_back(len_cct(&vec_tmp));
        err->push_back(err1);
        type->push_back(tp);
    }
    else if(LNRPSPLPRN(q_s,q_g,&vec_tmp,&t1,&u1,&v1,&err1,&tp))
    {
        vec_cct->push_back(vec_tmp);
        t->push_back(t1); u->push_back(u1); v->push_back(v1);
        length->push_back(len_cct(&vec_tmp));
        err->push_back(err1);
        type->push_back(tp);
    }
    else if(RPLNSNRNLP(q_s,q_g,&vec_tmp,&t1,&u1,&v1,&err1,&tp))
    {
        vec_cct->push_back(vec_tmp);
        t->push_back(t1); u->push_back(u1); v->push_back(v1);
        length->push_back(len_cct(&vec_tmp));
        err->push_back(err1);
        type->push_back(tp);
    }
    else if(RNLPSPRPLN(q_s,q_g,&vec_tmp,&t1,&u1,&v1,&err1,&tp))
    {
        vec_cct->push_back(vec_tmp);
        t->push_back(t1); u->push_back(u1); v->push_back(v1);
        length->push_back(len_cct(&vec_tmp));
        err->push_back(err1);
        type->push_back(tp);
    }
}
//---------------------------------------------------------------------------------------------------
//CCStateSpace::State CCStateSpace::CCPath::interpolate(double s, const CCStateSpace &space)
//s-step
void CCStateSpace::CCPath::vec_interpolate(double step, vector<CCStateSpace::State> *vec, const CCStateSpace &space)
{
    if(!cc_path.empty() && step >= CC_EPS)
    {
        for(size_t i=0;i<cc_path.size();++i)
            cc_path[i].vec_interpolate(step, vec, space);
    }
}

inline size_t min_index(vector<double> *d)
{
    size_t i=0;
    double t=(*d)[0];
    for(size_t j=1;j<d->size();++j)
    {   
        if(t>(*d)[j])
        {
            t=(*d)[j];
                i=j;
        }
    }
    return i;
}

CCStateSpace::CCPath::CCPath(CCStateSpace::State *start, CCStateSpace::State *goal, CCStateSpace *space):q_s(start),q_g(goal)
{
    vector<vector<CCStateSpace::CCTurn>> vec_cct;
    vector<double> tt; vector<double> uu; vector<double> vv;
    vector<double> l; vector<double> err1; vector<CCStateSpace::CCPathType> type1;
    t=-1;u=-1;v=-1;
    cc_length=0;
    cc_type=CC_NON;
    cc_err=-1;
    //
    space->TOPOLOGY(q_s, q_g, &vec_cct, &tt, &uu, &vv, &l, &err1, &type1);
    space->CSC(q_s, q_g, &vec_cct, &tt, &uu, &vv, &l, &err1, &type1);
    space->CCC(q_s, q_g, &vec_cct, &tt, &uu, &vv, &l, &err1, &type1);
    space->CCCC(q_s, q_g, &vec_cct, &tt, &uu, &vv, &l, &err1, &type1);
    space->CCSC(q_s, q_g, &vec_cct, &tt, &uu, &vv, &l, &err1, &type1);
    space->CCSCC(q_s, q_g, &vec_cct, &tt, &uu, &vv, &l, &err1, &type1);
    if(!vec_cct.empty())
    {
        size_t i=min_index(&l);
        t=tt[i];u=uu[i];v=vv[i];
        cc_length=l[i];
        cc_type=type1[i];
        cc_err=err1[i];
        cc_path=vec_cct[i];
    }
}

void CCStateSpace::interpolate(double *err, vector<CCStateSpace::State> *vec, CCStateSpace::State *q_s, CCStateSpace::State *q_g, int type, double *parameters, double step)
{

    const int *tp=CCStateSpace::PathType[type];
    vector<CCStateSpace::CCTurn> vec_cct;
    CCStateSpace::State *q=q_s;

    if(0<type && type<21)
    {
        for (int i=0; i<3; ++i)
        {
            if(tp[i]==1)
            {
                vec_cct.push_back(this->CC_LP(q, parameters[i]));
                q=&(vec_cct.back().q_g);
            }
            else if(tp[i]==2)
            {
                vec_cct.push_back(this->CC_RP(q, parameters[i]));
                q=&(vec_cct.back().q_g);
            }
            else if(tp[i]==3)
            {
                vec_cct.push_back(this->CC_LN(q, parameters[i]));
                q=&(vec_cct.back().q_g);
            }
            else if(tp[i]==4)
            {
                vec_cct.push_back(this->CC_RN(q, parameters[i]));
                q=&(vec_cct.back().q_g);
            }
            else if(tp[i]==5 && abs(parameters[i]>CC_EPS))
            {
                vec_cct.push_back(this->CC_SP(q, parameters[i]));
                q=&(vec_cct.back().q_g);
            }
            else if(tp[i]==6 && abs(parameters[i]>CC_EPS))
            {
                vec_cct.push_back(this->CC_SN(q, parameters[i]));
                q=&(vec_cct.back().q_g);
            }

        }
    }

    if(20<type && type<29)
    {
        double paras[4]={parameters[0],parameters[1],parameters[1],parameters[2]};
        for (int i=0; i<4; ++i)
        {
            if(tp[i]==1)
            {
                vec_cct.push_back(this->CC_LP(q, paras[i]));
                q=&(vec_cct.back().q_g);
            }
            else if(tp[i]==2)
            {
                vec_cct.push_back(this->CC_RP(q, paras[i]));
                q=&(vec_cct.back().q_g);
            }
            else if(tp[i]==3)
            {
                vec_cct.push_back(this->CC_LN(q, paras[i]));
                q=&(vec_cct.back().q_g);
            }
            else if(tp[i]==4)
            {
                vec_cct.push_back(this->CC_RN(q, paras[i]));
                q=&(vec_cct.back().q_g);
            }
            else if(tp[i]==5 && abs(parameters[i]>CC_EPS))
            {
                vec_cct.push_back(this->CC_SP(q, paras[i]));
                q=&(vec_cct.back().q_g);
            }
            else if(tp[i]==6 && abs(parameters[i]>CC_EPS))
            {
                vec_cct.push_back(this->CC_SN(q, paras[i]));
                q=&(vec_cct.back().q_g);
            }

        }
    }

    if((28<type && type<33) || (36<type && type<41))
    {
        double paras[4]={parameters[0],half_pi,parameters[1],parameters[2]};
        for (int i=0; i<4; ++i)
        {
            if(tp[i]==1)
            {
                vec_cct.push_back(this->CC_LP(q, paras[i]));
                q=&(vec_cct.back().q_g);
            }
            else if(tp[i]==2)
            {
                vec_cct.push_back(this->CC_RP(q, paras[i]));
                q=&(vec_cct.back().q_g);
            }
            else if(tp[i]==3)
            {
                vec_cct.push_back(this->CC_LN(q, paras[i]));
                q=&(vec_cct.back().q_g);
            }
            else if(tp[i]==4)
            {
                vec_cct.push_back(this->CC_RN(q, paras[i]));
                q=&(vec_cct.back().q_g);
            }
            else if(tp[i]==5 && abs(parameters[i]>CC_EPS))
            {
                vec_cct.push_back(this->CC_SP(q, paras[i]));
                q=&(vec_cct.back().q_g);
            }
            else if(tp[i]==6 && abs(parameters[i]>CC_EPS))
            {
                vec_cct.push_back(this->CC_SN(q, paras[i]));
                q=&(vec_cct.back().q_g);
            }
        }
    }

    if((32<type && type<37) || (40<type && type<45))
    {
        double paras[4]={parameters[0],parameters[1],half_pi,parameters[2]};
        for (int i=0; i<4; ++i)
        {
            if(tp[i]==1)
            {
                vec_cct.push_back(this->CC_LP(q, paras[i]));
                q=&(vec_cct.back().q_g);
            }
            else if(tp[i]==2)
            {
                vec_cct.push_back(this->CC_RP(q, paras[i]));
                q=&(vec_cct.back().q_g);
            }
            else if(tp[i]==3)
            {
                vec_cct.push_back(this->CC_LN(q, paras[i]));
                q=&(vec_cct.back().q_g);
            }
            else if(tp[i]==4)
            {
                vec_cct.push_back(this->CC_RN(q, paras[i]));
                q=&(vec_cct.back().q_g);
            }
            else if(tp[i]==5 && abs(parameters[i]>CC_EPS))
            {
                vec_cct.push_back(this->CC_SP(q, paras[i]));
                q=&(vec_cct.back().q_g);
            }
            else if(tp[i]==6 && abs(parameters[i]>CC_EPS))
            {
                vec_cct.push_back(this->CC_SN(q, paras[i]));
                q=&(vec_cct.back().q_g);
            }
        }
    }

    if(44<type && type<49)
    {
        double paras[5]={parameters[0],half_pi,parameters[1],half_pi,parameters[2]};
        for (int i=0; i<4; ++i)
        {
            if(tp[i]==1)
            {
                vec_cct.push_back(this->CC_LP(q, paras[i]));
                q=&(vec_cct.back().q_g);
            }
            else if(tp[i]==2)
            {
                vec_cct.push_back(this->CC_RP(q, paras[i]));
                q=&(vec_cct.back().q_g);
            }
            else if(tp[i]==3)
            {
                vec_cct.push_back(this->CC_LN(q, paras[i]));
                q=&(vec_cct.back().q_g);
            }
            else if(tp[i]==4)
            {
                vec_cct.push_back(this->CC_RN(q, paras[i]));
                q=&(vec_cct.back().q_g);
            }
            else if(tp[i]==5 && abs(parameters[i]>CC_EPS))
            {
                vec_cct.push_back(this->CC_SP(q, paras[i]));
                q=&(vec_cct.back().q_g);
            }
            else if(tp[i]==6 && abs(parameters[i]>CC_EPS))
            {
                vec_cct.push_back(this->CC_SN(q, paras[i]));
                q=&(vec_cct.back().q_g);
            }
        }
    }

    if(type == 49)
    {
        //1,two clothoids
        CCStateSpace::State q_s_ref(q_g->local_to_ref(q_s, *this));
        if(q_s_ref.theta > CC_EPS)
        {
            vec_cct.push_back(CCStateSpace::CCTurn(CCStateSpace::ccTurnType[10], *q_s, abs(q_s_ref.theta), parameters[0], *this)); //LN
        }
        else if(q_s_ref.theta < -CC_EPS)
        {
            vec_cct.push_back(CCStateSpace::CCTurn(CCStateSpace::ccTurnType[13], *q_s, abs(q_s_ref.theta), parameters[0], *this)); //RN
        }
        //2, line
        //CCStateSpace::State *q=q_s;
        CCStateSpace::State *q_r=&q_s_ref;
        if(!vec_cct.empty()) 
            {
                q=&(vec_cct.back().q_g);
                *q_r=q_g->local_to_ref(q, *this);
            }
        if(q_r->x < -CC_EPS)
        {
            vec_cct.push_back(CC_SP(q,-q_r->x));
        }
        else if(q_r->x > CC_EPS)
        {
            vec_cct.push_back(CC_SN(q, q_r->x));
        }
        //3, two clothoids + line + two clothoids
        if(!vec_cct.empty()) 
            {
                q=&(vec_cct.back().q_g);
                *q_r=q_g->local_to_ref(q, *this);
            }
        double l=abs(q_r->y)/tan(parameters[1]);
        if(q_r->y > CC_EPS)
        {
            //Lp
            vec_cct.push_back(CCStateSpace::CCTurn(CCStateSpace::ccTurnType[4], *q, 2*(parameters[1]), parameters[2], *this));
            //Sn
            vec_cct.push_back(CCStateSpace::CCTurn(CCStateSpace::ccTurnType[1], (vec_cct.back()).q_g, l,*this));
            //Rp
            vec_cct.push_back(CCStateSpace::CCTurn(CCStateSpace::ccTurnType[7], (vec_cct.back()).q_g, 2*(parameters[1]), parameters[2], *this));
        }
        else if(q_r->y < -CC_EPS)
        {
            //Rp
            vec_cct.push_back(CCStateSpace::CCTurn(CCStateSpace::ccTurnType[7], *q, 2*(parameters[1]), parameters[2], *this));
            //Sn
            vec_cct.push_back(CCStateSpace::CCTurn(CCStateSpace::ccTurnType[1], (vec_cct.back()).q_g, l,*this));
            //Lp
            vec_cct.push_back(CCStateSpace::CCTurn(CCStateSpace::ccTurnType[4], (vec_cct.back()).q_g, 2*(parameters[1]), parameters[2], *this));
        }
    }

    if(!vec_cct.empty() && step >= CC_EPS)
    {
        for(size_t i=0;i<vec_cct.size();++i)
            vec_cct[i].vec_interpolate(step, vec, *this);
    }
    if(!vec_cct.empty()) 
        q=&(vec_cct.back().q_g);
    *err=this->distance(q, q_g);
}
