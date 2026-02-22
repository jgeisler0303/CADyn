#ifndef RK1CONDENSED_HPP_
#define RK1CONDENSED_HPP_

#include <exception>
#include <string>
#include <iosfwd>
#include <cmath>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include "OptionInfo.hpp"

class RK1condensedException: public std::exception {
public:
    RK1condensedException(const std::string& msg= "General RK1condensed integration exception") :
        m_msg(msg)
    {  }
    
    virtual const char* what() const throw () {
        return m_msg.c_str();
    }
    
    const std::string m_msg;
};

template<typename Derived>
inline bool is_finite(const Eigen::MatrixBase<Derived>& x)
{
	return ( (x - x).array() == (x - x).array()).all();
}

// This is an example stub definition of what a problem class needs to implement to be used with the condensed implicit RK1 solver
/*
class ProblemDefinition {
public:
    static constexpr int dimI = 3;              // number of states that will be condensed out
    static constexpr int dimII = 4;             // number of states that will be kept
    static constexpr int dimX = dimI + dimII;   // number of states
    static constexpr int dimIn = 1;             // number of inputs
    static constexpr int dimOut = 2;            // number of outputs

    typedef Eigen::Matrix<real_type, dimI+dimII, 1> VecX;
    typedef Eigen::Matrix<real_type, dimI+dimII, dimI+dimII> MatX;
    typedef Eigen::Matrix<real_type, dimII, 1> VecII;
    typedef Eigen::Matrix<real_type, dimI, dimII> MatI_II;
    typedef Eigen::Matrix<real_type, dimII, dimI> MatII_I;
    typedef Eigen::Matrix<real_type, dimII, dimII> MatII;
    typedef Eigen::Matrix<real_type, dimIn, 1> VecIn;
    typedef Eigen::Matrix<real_type, dimI+dimII, dimIn> MatIn;
    typedef Eigen::Matrix<real_type, dimII, dimIn> MatII_In;
    typedef Eigen::Matrix<real_type, dimOut, 1> VecOut;
    typedef Eigen::Matrix<real_type, dimOut, dimI+dimII> MatOutX;
    typedef Eigen::Matrix<real_type, dimOut, dimII> MatOutII;
    typedef Eigen::Matrix<real_type, dimOut, dimIn> MatOutIn;

    const MatI_II E = Eigen::Map<const MatI_II>(
        (real_type[]){  1.0,0.0,0.0,
                        0.0,1.0,0.0});      

    // residual evaluation
    virtual void evaluateResidual(
        const VecX& x,
        const VecII& xdot_II,
        const VecIn& u,
        real_type t) = 0;

    // Jacobian blocks for f_II written to class members df_II_dx_I, df_II_dx_II, df_II_dxdot_II
    virtual void evaluateResidualAndJacobians() = 0;
    virtual void evaluateInputJacobian() = 0;
    virtual void evaluateOutput() = 0;
    virtual void evaluateOutputJacobian() = 0;
    virtual void precalcConsts();

    typename inputs_t inputs;
    typename outputs_t outputs;
    typename states_t states;

    // Inputs and results
    real_type t = 0.0;  // current time
    VecX x;             // state at evaluation time
    VecII xdot_II;      // state derivative at evaluation time
    VecIn u;
    VecOut y;

    // evaluation results
    VecII f_II;
    MatII_I df_II_dx_I;
    MatII df_II_dx_II;
    MatII df_II_dxdot_II;
    MatII_In df_II_du;

    VecOut y;
    MatOutX dy_dx;
    MatOutII dy_ddotx_II;
    MatOutIn dy_du;
    
    Parameters param;
};
*/

template <class problem_class>
class RK1condensed : public OptionsAccessor, public problem_class {
public:
    static const int nI  = problem_class::dimI;
    static const int nII = problem_class::dimII;
    static const int nX = nI + nII;
    static const int nIn = problem_class::dimIn;
    static const int nOut = problem_class::dimOut;

    // Access typedefs from template parameter
    using real_type = typename problem_class::real_type;
    using VecI = typename problem_class::VecI;
    using VecII = typename problem_class::VecII;
    using VecX = typename problem_class::VecX;
    using MatI = typename problem_class::MatI;
    using MatII = typename problem_class::MatII;
    using MatII_I = typename problem_class::MatII_I;
    using MatX = typename problem_class::MatX;
    using MatIn = typename problem_class::MatIn;
    using VecIn = typename problem_class::VecIn;
    using MatII_In = typename problem_class::MatII_In;
    using VecOut = typename problem_class::VecOut;
    using MatOutX = typename problem_class::MatOutX;
    using MatOutII = typename problem_class::MatOutII;
    using MatOutIn = typename problem_class::MatOutIn;

    // Make problem class members available in this class
    using problem_class::name;
    using problem_class::description;

    using problem_class::inputs;
    using problem_class::outputs;
    using problem_class::states;

    using problem_class::in_name;
    using problem_class::out_name;
    using problem_class::state_name;

    using problem_class::inputs_idx;
    using problem_class::outputs_idx;
    using problem_class::states_idx;

    using problem_class::param;
    using problem_class::t;
    using problem_class::x;
    using problem_class::xdot;
    using problem_class::u;
    using problem_class::y;
    using problem_class::f_II;
    using problem_class::E;
    using problem_class::df_II_dx_I;
    using problem_class::df_II_dx_II;
    using problem_class::df_II_dxdot_II;
    using problem_class::df_II_du;
    using problem_class::dy_dx;
    using problem_class::dy_dxdot;
    using problem_class::dy_du;

    // Make problem class methods available in this class
    using problem_class::evaluateResidual;
    using problem_class::evaluateResidualAndJacobians;
    using problem_class::evaluateInputJacobian;
    using problem_class::evaluateOutput;
    using problem_class::evaluateOutputJacobian;
    using problem_class::precalcConsts;

    RK1condensed() {
        locked_states.fill(false);

        options_map["AbsTol"]= OptionInfo(&AbsTol, 1e-12, std::numeric_limits<real_type>::infinity());
        options_map["RelTol"]= OptionInfo(&RelTol, 1e-12, std::numeric_limits<real_type>::infinity());
        options_map["StepTol"]= OptionInfo(&StepTol, 1e-12, std::numeric_limits<real_type>::infinity());
        options_map["StepAdjMin"]= OptionInfo(&StepAdjMin, 0.0, 1.0);
        options_map["StepAdjMax"]= OptionInfo(&StepAdjMax, 1.0, std::numeric_limits<real_type>::infinity());
        options_map["StepAdjSafe"]= OptionInfo(&StepAdjSafe, 0.0, 1.0);
        options_map["hminmin"]= OptionInfo(&hminmin, 1e-12, std::numeric_limits<real_type>::infinity());
        options_map["jac_recalc_step"]= OptionInfo(&jac_recalc_step, 1.0, std::numeric_limits<real_type>::infinity());
        options_map["max_steps"]= OptionInfo(&max_steps, 1.0, std::numeric_limits<real_type>::infinity());
        options_map["rk_a"]= OptionInfo(&rk_a, 0.0, 1.0);
    }

    enum step_result {
        success,
        max_iterations_reached,
        err_not_finite
    };

    step_result oneStep(real_type h, bool hmodified= true, bool add_damping= false);
    real_type computeError(real_type h);
    void interval(real_type tfinal, real_type &h_save, real_type hmax, bool pseudo_time_continuation= false);
    void intervalWithSens(real_type ts) { return intervalWithSens(ts, ts); };
    void intervalWithSens(real_type ts, real_type h)  {
        interval(t+ts, h, ts);
        sensitivities(ts);
    };
    void sensitivities(real_type ts);

    void staticEquilibrium();
    void staticEquilibriumWithLin();

    MatX fx;
    MatIn fu;
    MatOutX gx;
    MatOutIn gu;
    
    // Solver data
    Eigen::FullPivLU<MatII> LU;
    VecX x_last;        // last state before step (x_k) for error estimation

    int iter;
    int sub_iter;
    int n_back_steps;
    int n_steps;

    // Options
    real_type AbsTol= 1e-6;
    real_type RelTol= 1e-6;
    real_type StepTol= 1e-8;
    real_type StepAdjMin= 0.1;
    real_type StepAdjMax= 5.0;
    real_type StepAdjSafe= 0.9;
    real_type hminmin= 1E-8;
    int jac_recalc_step= 4;
    int max_steps= 10;
    real_type rk_a= 1.0; // default: implicit Euler

    std::array<bool, nX> locked_states;
};

template <class problem_class>
typename RK1condensed<problem_class>::step_result RK1condensed<problem_class>::oneStep(real_type h, bool hmodified, bool add_damping) {
    real_type ha = h*rk_a;
    VecX x_store = x;
    VecII xdot_store = xdot;
    real_type t_store = t;
    t += ha;

    int iter_jac= 0;
    if(hmodified && jac_recalc_step>1) iter_jac= 1;

    real_type err;
    VecII Dii = VecII::Zero();
    for (iter = 1; iter <= max_steps; ++iter, ++sub_iter) {
        // x_{k+1} = x_k + a*h*k
        // for x_II, locked states are already considered because xdot is zero for locked states
        x.tail(nII) = x_store.tail(nII) + ha * xdot;
        // compute x_II first because k_II also affects x_I
        for (int i = 0; i < nI; ++i) {
            if(locked_states[i])
                x[i] = x_store[i];
            else
                x[i] = x_store[i]  + ha * E.row(i) * x.tail(nII);
        }


        if((iter%jac_recalc_step)==iter_jac) {
            evaluateResidualAndJacobians();

            // Construct reduced blocks
            // A = a*h*df_II_dx_I,  df_II_dxdot_I = 0
            // B = a*h*df_II_dx_II + df_II_dxdot_II
            // (B + A a*h E) Δk_II = -f_II
            MatII Jacobian = df_II_dxdot_II + ha*df_II_dx_II + ha*ha*df_II_dx_I*E;
            if(add_damping) {
                VecII Keffdiag= (df_II_dx_I*E).diagonal();
                for (int i = 0; i < nII; ++i) {
                    real_type m = std::max(df_II_dxdot_II(i,i), 1e-12);
                    real_type k = std::max(std::abs(Keffdiag(i)), 1e-12);
                    Dii[i] = 2.0 * std::sqrt(m * k);
                }
                Jacobian += ha*Dii.asDiagonal();
            }
            for (int i = 0; i < nII; ++i) {    
                if (locked_states[nI + i]) {
                    Jacobian.col(i).setZero();
                    Jacobian.row(i).setZero();
                    Jacobian(i, i) = 1.0;
                }
            }
            LU.compute(Jacobian);
        } else {
            evaluateResidual();
        }

        VecII residual = f_II;
        if(add_damping) {
            residual += Dii.cwiseProduct(x.tail(nII));
        }
        for (int i = 0; i < nII; ++i) {
            if(locked_states[nI + i]) {
                residual[i] = 0.0;
            }
        }
        VecII Dxdot = LU.solve(-residual);

        if(!is_finite(Dxdot))
            throw RK1condensedException("Step calculation: Newton step is not finite.");

        xdot += Dxdot;

        err= Dxdot.norm() / (sqrt(1.0*nII) * (1.0 + xdot.norm()));
        if (err < StepTol)
            break;
    }

    if(iter >= max_steps || !std::isfinite(err)) {
            // restore state to before the failed step
            t= t_store;
            x= x_store;
            xdot= xdot_store;
            if(iter >= max_steps)
                return max_iterations_reached;
            else
                return err_not_finite;
    }

    // take the full step
    t = t_store + h;
    x_last = x_store; // save the state before the step (x_k) for error estimation
    VecII x_II_ha_xdot = x_store.tail(nII) + ha * xdot;
    for (int i = 0; i < nI; ++i) {
        if(locked_states[i])
            x[i] = x_store[i];
        else
            x[i] = x_store[i]  + h * E.row(i) * x_II_ha_xdot;
    }
    x.tail(nII) = x_store.tail(nII) + h * xdot;

    return success;
}

template <class problem_class>
typename RK1condensed<problem_class>::real_type RK1condensed<problem_class>::computeError(real_type h) {
    // If rk_a == 0.5, the method is already second order and the error estimate will be zero. So this error estimate is not meaningful and we should throw an error 
    // The neglected higher order terms in the Taylor expansion will also be relevant for rk_a around 0.5.
    if(std::abs(rk_a - 0.5) < 0.1) {
        throw RK1condensedException("Error estimation is not meaningful for rk_a close to 0.5. Please set rk_a to a value sufficiently different from 0.5 (e.g. 0.25 or 0.75) to get a meaningful error estimate.");
    }
    // The error estimate is based on the local truncation error of the method
    VecX tau = VecX::Zero();

    // We first consider a non partitioned, and explicit system for simplicity:
    //   dot x = g(x, t)
    // The RK1 a-method is: x(t+h) = x(t) + h k
    //   where k = g(x + a h k, t + a h)
    //   if g is expanded about x_k, t_k: g approx g + dg/dx Delta x + dg/dt Delta t
    //   with Delta x = a h k and Delta t = a h, we have: 
    //     k approx g + a h (dg/dx k + dg/dt)
    //   substituting k = g and then dg/dx g + dg/dt = dot g (= dg/dx dotx + dg/dt) yields
    //     k = g + a h dot g + O(h^2)
    // The RK1 method expanded is thus: x(t+h) = x(t) + h g + h^2 a dot g + O(h^3)
    // The exact solution expanded is:  x(t+h) = x(t) + h g + h^2/2 dot ​g​ + O(h3)
    // And we arrive at the local truncation error estimate tau approx h^2 (1/2-a) dot k

    // Estimate local truncation error for part II
    // To get dot k we start With k from the implicit equation: f(x, k, t) = 0 and
    //   differentiate w.r.t. time: df/dx * xdot + df/dxdot * dot k + df/dt = 0
    //   thus dot k = -(df/dxdot)^-1 * (df/dx * xdot + df/dt)
    // For the partitioned system we have: f_II(x_I, x_II, xdot_II, u, t) = 0,
    //   and with dot x_I = E xdot_II,
    // we thus get: dot k_II approx -(df_II_dxdot_II)^-1 * (df_II_dx_I * E * xdot_II + df_II_dx_II * xdot_II + df_II_dt)
    VecII x_II_stage = x_last.tail(nII) + rk_a*h*xdot;
    VecII rhs = df_II_dx_I * E * x_II_stage + df_II_dx_II * xdot; // + df_II_dt;

    // To avoid inverting df_II_dxdot_II, we can use the Jacobian and its LU from the Newton iteration (especiall if h is small)
    // VecII kdot = LU.solve(-rhs);
    //   or we can use a diagonal approximation of df_II_dxdot_II for
    VecII kdot = -rhs.cwiseQuotient(df_II_dxdot_II.diagonal().cwiseMax(1e-12));
    
    // assignment in loop to avoid AVX -Warray-bounds warning for 3-element vectors.
    const real_type err_fac = h*h*(0.5 - rk_a);
    for (int i = 0; i < nII; ++i) {
        tau[nI + i] = err_fac * kdot[i];
    }

    // Estimate local truncation error for condensed part (I)
    // Taylor expansion about t_k: x_I​(t_k​+h) = x_{I,k​} + h E x_{II,k​} + h^2/2 ​E dot x_{II,k}​ + O(h3)
    // Method of this solver:      x_I​(t_k​+h) = x_{I,k}​ + h E x_{II,k} ​+ h^2 a ​E k_II
    // Subtract numerical from exact: tau_I = h^2 E (1/2 dot x_{II,k}​ - a k_II) + O(h^3)
    // Substitute dot x_{II,k}​= k_II + O(h): tau_I approx h^2/2 E (1- a/2) dot x_{II,k} = E tau_II 
    VecI tau_I = err_fac * E * xdot;
    for (int i = 0; i < nI; ++i) {
        tau[i] = tau_I[i];
    }

    // reconstruct full xdot for scaling the error
    VecX xdot_full;
    xdot_full.tail(nII) = xdot;
    xdot_full.head(nI)  = E * x.tail(nII);

    real_type sum = 0.0;

    for (int i = 0; i < nX; ++i) {
        real_type xi = std::abs(x[i]);
        real_type xdoti = std::abs(xdot_full[i]);

        // Dont' take max of AbsTol and RelTol*|x| to make behaviour smoother. This is standard practice in ODE solvers.
        real_type sc = AbsTol + RelTol * std::max(xi, h * xdoti);
        
        if(!locked_states[i]) {
            real_type r = tau[i] / sc;
            sum += r * r;
        }
    }

    return std::sqrt(sum / nX);
}

template <class problem_class>
void RK1condensed<problem_class>::interval(real_type tfinal, real_type &h_save, real_type hmax, bool pseudo_time_continuation) {
    // This function implements a simple adaptive time stepping loop around oneStep() with error control based on computeError()
    // and step size adaptation similar to standard ODE solvers. It also keeps track of the number of steps taken, number of back steps, etc. for performance monitoring.
    // Integrate over one interval until tfinal. Start with step size h_save. When longer steps are possible, never increase step size greater than hmax.
    bool hchanged= true;
    n_steps= 0;         // number of successful steps taken (calls to oneStep)
    sub_iter= 0;        // number of calls to oneStep (including back steps
    n_back_steps= 0;    // number of steps that had to be redone with smaller step size

    if(!is_finite(xdot) || !is_finite(x) || !is_finite(u))
        throw RK1condensedException("Interval calculation: initial values of x, xdot, u are not finite.");
    
    if (h_save>hmax) h_save= hmax;
    while(t < tfinal) {
        // Adjust step size to hit tfinal exactly
        if((t+h_save) > tfinal) {
            h_save= tfinal-t;
            hchanged= true;
        } else if((t+2.1*h_save) > tfinal) {
            // If the next two steps would take us past tfinal, we need to adjust the step size to avoid one final very small step at the end. This is a common practice in ODE solvers to improve performance and avoid numerical issues with very small steps.
            h_save= 0.5*(tfinal-t);
            hchanged= true;
        }
        // store current state in case we need to redo the step with smaller step size
        real_type time_store= t;
        VecX x_store= x;
        VecII xdot_store= xdot;

        RK1condensed<problem_class>::step_result res = oneStep(h_save, hchanged, pseudo_time_continuation);
        real_type err= computeError(h_save);
        n_steps++;

        // If step failed or error is too large, reduce step size and redo the step
        if((res == max_iterations_reached) || (res == err_not_finite) || (err>1.0))   {
            // Check if step size can be adjusted based on error estimate or if we need to reduce it more drastically due to failure of oneStep
            if((res == max_iterations_reached) || (res == err_not_finite) || !std::isfinite(err))
                h_save*= 0.25;
            else
                // error is already scaled by tolernces in computeError, so we can just use sqrt(1/err) for step size adaptation
                h_save*= std::max(StepAdjSafe*sqrt(1.0/err), StepAdjMin);
            hchanged= true;
            
            // restore state to before the failed step
            t= time_store;
            x= x_store;
            xdot= xdot_store;
            
            if(h_save<hminmin)
                throw RK1condensedException("Interval calculation: dynamic step size too small.");
            
            n_back_steps++;
        } else if (err < StepAdjSafe*StepAdjSafe) {
            h_save*= std::min(StepAdjSafe*sqrt(1.0/err), StepAdjMax);
            if(h_save>hmax) h_save= hmax;
            hchanged= true;
        } else {
            hchanged= false;
        }

        if(pseudo_time_continuation) {
            if(n_steps>10000)
                throw RK1condensedException("Interval calculation: no static equilibrium found after 10000 steps.");

            real_type stationarity_check = 0.0;

            VecI dot_xI = E*x.tail(nII);
            for (int i = 0; i < nI; ++i) {
                if (!locked_states[i]) {
                    stationarity_check += dot_xI[i] * dot_xI[i];
                }
            }
            for (int i = 0; i < nII; ++i) {
                if (!locked_states[nI + i]) {
                    stationarity_check += xdot[i] * xdot[i];
                }
            }
            if(err < 1.0 && std::sqrt(stationarity_check) < StepTol)
                break;
        }
    }
}

template <class problem_class>
void RK1condensed<problem_class>::sensitivities(real_type h) {
    // compute sensitivities using the implicit function theorem. We have f_II(x, k_II, u)= 0 at the solution point. We want to compute dy/du.
    evaluateResidualAndJacobians();
    evaluateInputJacobian();
    evaluateOutputJacobian();
    
    real_type ha= h*rk_a;
    
    // partial f_II,ev / partial k_II: reused Jacobian from interval/oneStep calculation
    // partial k_II / partial x_I,k
    MatII_I Pk_II_Px_Ik= -LU.solve(df_II_dx_I);
    // partial k_II / partial x_II,k
    MatII Pk_II_Px_IIk= -LU.solve(df_II_dx_II + ha*df_II_dx_I*E);
    // partial k_II / partial u_k
    MatII_In Pk_II_Pu= -LU.solve(df_II_du);
    
    // partial x_I,k+1 / partial x_I,k = I + ha*h*E*partial k_II / partial x_I,k
    fx.block(0, 0, nI, nI)= MatI::Identity() + ha*h*E*Pk_II_Px_Ik;
    // partial x_I,k+1 / partial x_II,k = h*E + ha*h*E*partial k_II / partial x_II,k
    fx.block(0, nI, nI, nII)= h*E + ha*h*E*Pk_II_Px_IIk;
    // partial x_II,k+1 / partial x_I,k = h*partial k_II / partial x_I,k
    fx.block(nI, 0, nII, nI)= h*Pk_II_Px_Ik;
    // partial x_II,k+1 / partial x_II,k = I + h*partial k_II / partial x_II,k
    fx.block(nI, nI, nII, nII)= MatII::Identity() + h*Pk_II_Px_IIk;

    // partial x_I,k+1 / partial u_k = ha*h*E*partial k_II / partial u_k
    fu.block(0, 0, nI, nIn)= ha*h*E*Pk_II_Pu;
    // partial x_II,k+1 / partial u_k = h*partial k_II / partial u_k
    fu.block(nI, 0, nII, nIn)= h*Pk_II_Pu;
    
    gx = dy_dx;
    gx.block(0, 0, nOut, nI)      += dy_dxdot*Pk_II_Px_Ik;
    gx.block(0, nI, nOut, nII)    += dy_dxdot*Pk_II_Px_IIk;
    gu = dy_du  + dy_dxdot*Pk_II_Pu;
}

template <class problem_class>
void RK1condensed<problem_class>::staticEquilibrium() {
    // use pseudo time continuation with interval() with added damping to find a static equilibrium point.
    real_type rk_a_save= rk_a;
    rk_a= 1.0; // use fully implicit Euler for better convergence of static equilibrium calculation
    // take one large step and terminate when k_II is almost zero
    real_type h_eq = 1.0;
    interval(std::numeric_limits<real_type>::infinity(), h_eq, std::numeric_limits<real_type>::infinity(), true);
    rk_a= rk_a_save;
    t = 0.0; // reset time to zero after pseudo time continuation
}

template <class problem_class>
void RK1condensed<problem_class>::staticEquilibriumWithLin() {
    staticEquilibrium();
    evaluateResidualAndJacobians();
    evaluateInputJacobian();
    evaluateOutputJacobian();
    
    // return continuous time linearization
    fx.block(0, 0, nI, nI).setZero();
    fx.block(0, nI, nI, nII)= E;

    LU.compute(df_II_dxdot_II);
    fx.block(nI, 0, nII, nI)= -LU.solve(df_II_dx_I);
    fx.block(nI, nI, nII, nII)= -LU.solve(df_II_dx_II);
    
    fu.block(0, 0, nI, nIn).setZero();
    fu.block(nI, 0, nII, nIn)= -LU.solve(df_II_du);
}    

#endif /* RK1CONDENSED_HPP_ */
