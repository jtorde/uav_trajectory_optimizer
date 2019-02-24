/****************************************************************************
 *   Copyright (c) 2019 Jesus Tordesillas Torres. All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in
 *    the documentation and/or other materials provided with the
 *    distribution.
 * 3. Neither the name snap nor the names of its contributors may be
 *    used to endorse or promote products derived from this software
 *    without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
 * FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 * COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
 * INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
 * BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS
 * OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
 * AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 * LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
 * ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 *
 ****************************************************************************/
#ifndef SOLVERS_HPP
#define SOLVERS_HPP
#include <Eigen/Dense>
#include <iostream>

#include "cvxgen/interface_vel.h"
#include "cvxgen/interface_accel.h"
#include "cvxgen/interface_jerk.h"

#include <unsupported/Eigen/Polynomials>

#define VEL 1
#define ACCEL 2
#define JERK 3

#define N_VEL 15
#define N_ACCEL 15
#define N_JERK 10

template <int INPUT_ORDER>
class Solver
{
public:
  Solver();
  void set_x0(double x0[]);
  void set_xf(double xf[]);
  void set_max(double max_values[INPUT_ORDER]);

  void genNewTraj();

  int getN();
  void setq(double q);
  double** getState();
  double** getInput();
  double getCost();

protected:
  bool checkConvergence(double xf_opt[3 * INPUT_ORDER]);
  void callOptimizer();
  double getDTInitial();
  double dt_;  // time step found by the solver
  int N_;
  double xf_[3 * INPUT_ORDER];
  double x0_[3 * INPUT_ORDER];
  double v_max_;
  double a_max_;
  double j_max_;
  double q_;  // weight to the 2nd term in the cost function
  double** x_;
  double** u_;
};

template <int INPUT_ORDER>
double Solver<INPUT_ORDER>::getCost()
{
  double cost;
  switch (INPUT_ORDER)
  {
    case VEL:
      cost = vel_get_cost();
      break;
    case ACCEL:
      cost = accel_get_cost();
      break;
    case JERK:
      cost = jerk_get_cost();
      break;
  }
}

template <int INPUT_ORDER>
double** Solver<INPUT_ORDER>::getState()
{
  return x_;
}

template <int INPUT_ORDER>
double** Solver<INPUT_ORDER>::getInput()
{
  return u_;
}

template <int INPUT_ORDER>
Solver<INPUT_ORDER>::Solver()
{
  v_max_ = 20;
  a_max_ = 2;
  j_max_ = 20;
  switch (INPUT_ORDER)
  {
    case VEL:
      vel_initialize_optimizer();
      N_ = N_VEL;
      break;
    case ACCEL:
      accel_initialize_optimizer();
      N_ = N_ACCEL;
      break;
    case JERK:
      jerk_initialize_optimizer();
      N_ = N_JERK;
      break;
  }
}

template <int INPUT_ORDER>
int Solver<INPUT_ORDER>::getN()
{
  return N_;
}

template <int INPUT_ORDER>
void Solver<INPUT_ORDER>::setq(double q)
{
  q_ = q;
}

template <int INPUT_ORDER>
void Solver<INPUT_ORDER>::set_x0(double x0[])
{
  for (int i = 0; i < 3 * INPUT_ORDER; i++)
  {
    x0_[i] = x0[i];
  }
}

template <int INPUT_ORDER>
void Solver<INPUT_ORDER>::set_xf(double xf[])
{
  for (int i = 0; i < 3 * INPUT_ORDER; i++)
  {
    xf_[i] = xf[i];
  }
}

template <int INPUT_ORDER>
void Solver<INPUT_ORDER>::set_max(double max_values[INPUT_ORDER])
{
  switch (INPUT_ORDER)
  {
    case VEL:
      v_max_ = max_values[0];
      break;
    case ACCEL:
      v_max_ = max_values[0];
      a_max_ = max_values[1];
      break;
    case JERK:
      v_max_ = max_values[0];
      a_max_ = max_values[1];
      j_max_ = max_values[2];
      break;
  }
}

template <int INPUT_ORDER>
bool Solver<INPUT_ORDER>::checkConvergence(double xf_opt[3 * INPUT_ORDER])
{
  bool converged = false;
  float d2 = 0;   // distance in position squared
  float dv2 = 0;  // distance in velocity squared
  float da2 = 0;  // distance in acceleration squared

  switch (INPUT_ORDER)
  {
    case VEL:
      for (int i = 0; i < 3; i++)
      {
        d2 += pow(xf_[i] - xf_opt[i], 2);
      }
      converged = (sqrt(d2) < 0.2) ? true : false;
      break;
    case ACCEL:
      for (int i = 0; i < 3; i++)
      {
        d2 += pow(xf_[i] - xf_opt[i], 2);
        dv2 += pow(xf_[i + 3] - xf_opt[i + 3], 2);
      }
      converged = (sqrt(d2) < 0.2 && sqrt(dv2) < 0.2) ? true : false;
      break;
    case JERK:
      for (int i = 0; i < 3; i++)
      {
        d2 += pow(xf_[i] - xf_opt[i], 2);
        dv2 += pow(xf_[i + 3] - xf_opt[i + 3], 2);
        da2 += pow(xf_[i + 6] - xf_opt[i + 6], 2);
      }
      converged = (sqrt(d2) < 0.2 && sqrt(dv2) < 0.2 && 1) ? true : false;
      break;
  }

  return converged;
}

template <int INPUT_ORDER>
void Solver<INPUT_ORDER>::genNewTraj()
{
  callOptimizer();

  switch (INPUT_ORDER)
  {
    case VEL:
      x_ = vel_get_state();
      u_ = vel_get_control();

      break;
    case ACCEL:
      x_ = accel_get_state();
      u_ = accel_get_control();

      break;
    case JERK:
      x_ = jerk_get_state();
      u_ = jerk_get_control();

      break;
  }
}

template <int INPUT_ORDER>
void Solver<INPUT_ORDER>::callOptimizer()
{
  bool converged = false;

  double dt = getDTInitial();

  double** x;
  int i = 0;
  int r = 0;
  while (1)
  {
    dt = dt + 3 * 0.05;
    i = i + 1;
    switch (INPUT_ORDER)
    {
      case VEL:
      {
        vel_load_default_data(dt, v_max_, x0_, xf_, q_);
        r = vel_optimize();
        if (r == 1)
        {
          x = vel_get_state();
          converged = checkConvergence(x[N_]);
        }
        break;
      }
      case ACCEL:
      {
        accel_load_default_data(dt, v_max_, a_max_, x0_, xf_, q_);
        r = accel_optimize();
        if (r == 1)
        {
          x = accel_get_state();
          converged = checkConvergence(x[N_]);
        }
        break;
      }
      case JERK:
      {
        jerk_load_default_data(dt, v_max_, a_max_, j_max_, x0_, xf_, q_);
        r = jerk_optimize();
        if (r == 1)
        {
          x = jerk_get_state();
          converged = checkConvergence(x[N_]);
        }
        break;
      }
    }
    if (converged == 1)
    {
      break;
    }
  }

  if (i > 1)
  {
    printf("Iterations = %d\n", i);
    printf("Iterations>1, if you increase dt at the beginning, it would be faster\n");
  }
  dt_ = dt;
}

inline double MinPositiveElement(std::vector<double> v)
{
  std::sort(v.begin(), v.end());  // sorted in ascending order
  double min_value = 0;
  for (int i = 0; i < v.size(); i++)
  {
    if (v[i] > 0)
    {
      min_value = v[i];
      break;
    }
  }
  return min_value;
}

template <int INPUT_ORDER>
double Solver<INPUT_ORDER>::getDTInitial()
{
  double dt_initial = 0;
  float t_vx = 0;
  float t_vy = 0;
  float t_vz = 0;
  float t_ax = 0;
  float t_ay = 0;
  float t_az = 0;
  float t_jx = 0;
  float t_jy = 0;
  float t_jz = 0;

  t_vx = (xf_[0] - x0_[0]) / v_max_;
  t_vy = (xf_[1] - x0_[1]) / v_max_;
  t_vz = (xf_[2] - x0_[2]) / v_max_;

  switch (INPUT_ORDER)
  {
    case JERK:
    {
      float jerkx = copysign(1, xf_[0] - x0_[0]) * j_max_;
      float jerky = copysign(1, xf_[1] - x0_[1]) * j_max_;
      float jerkz = copysign(1, xf_[2] - x0_[2]) * j_max_;
      float a0x = x0_[6];
      float a0y = x0_[7];
      float a0z = x0_[8];
      float v0x = x0_[3];
      float v0y = x0_[4];
      float v0z = x0_[5];

      // polynomial ax3+bx2+cx+d=0 --> coeff=[d c b a]
      Eigen::Vector4d coeff_jx(x0_[0] - xf_[0], v0x, a0x / 2.0, jerkx / 6.0);
      Eigen::Vector4d coeff_jy(x0_[1] - xf_[1], v0y, a0y / 2.0, jerky / 6.0);
      Eigen::Vector4d coeff_jz(x0_[2] - xf_[2], v0z, a0z / 2.0, jerkz / 6.0);

      /*  std::cout << "Coefficients for jerk" << std::endl;
        std::cout << "Coeffx=" << coeff_jx.transpose() << std::endl;
        std::cout << "Coeffy=" << coeff_jy.transpose() << std::endl;
        std::cout << "Coeffz=" << coeff_jz.transpose() << std::endl;*/

      Eigen::PolynomialSolver<double, Eigen::Dynamic> psolve_jx(coeff_jx);
      Eigen::PolynomialSolver<double, Eigen::Dynamic> psolve_jy(coeff_jy);
      Eigen::PolynomialSolver<double, Eigen::Dynamic> psolve_jz(coeff_jz);

      std::vector<double> realRoots_jx;
      std::vector<double> realRoots_jy;
      std::vector<double> realRoots_jz;
      psolve_jx.realRoots(realRoots_jx);
      psolve_jy.realRoots(realRoots_jy);
      psolve_jz.realRoots(realRoots_jz);

      t_jx = MinPositiveElement(realRoots_jx);
      t_jy = MinPositiveElement(realRoots_jy);
      t_jz = MinPositiveElement(realRoots_jz);

      // printf("Times: t_jx, t_jy, t_jz:\n");
      // std::cout << t_jx << "  " << t_jy << "  " << t_jz << std::endl;

      // Here there is no a break
    }
    case ACCEL:
    {
      float accelx = copysign(1, xf_[0] - x0_[0]) * a_max_;
      float accely = copysign(1, xf_[1] - x0_[1]) * a_max_;
      float accelz = copysign(1, xf_[2] - x0_[2]) * a_max_;
      float v0x = x0_[3];
      float v0y = x0_[4];
      float v0z = x0_[5];

      // polynomial ax2+bx+c=0 --> coeff=[c b a]
      Eigen::Vector3d coeff_ax(x0_[0] - xf_[0], v0x, 0.5 * accelx);
      Eigen::Vector3d coeff_ay(x0_[1] - xf_[1], v0y, 0.5 * accely);
      Eigen::Vector3d coeff_az(x0_[2] - xf_[2], v0z, 0.5 * accelz);

      Eigen::PolynomialSolver<double, Eigen::Dynamic> psolve_ax(coeff_ax);
      Eigen::PolynomialSolver<double, Eigen::Dynamic> psolve_ay(coeff_ay);
      Eigen::PolynomialSolver<double, Eigen::Dynamic> psolve_az(coeff_az);

      std::vector<double> realRoots_ax;
      std::vector<double> realRoots_ay;
      std::vector<double> realRoots_az;
      psolve_ax.realRoots(realRoots_ax);
      psolve_ay.realRoots(realRoots_ay);
      psolve_az.realRoots(realRoots_az);

      t_ax = MinPositiveElement(realRoots_ax);
      t_ay = MinPositiveElement(realRoots_ay);
      t_az = MinPositiveElement(realRoots_az);
    }
    case VEL:
    {
      // I'm done
      break;
    }
  }
  dt_initial = std::max({ t_vx, t_vy, t_vz, t_ax, t_ay, t_az, t_jx, t_jy, t_jz }) / N_;
  if (dt_initial > 10000)  // happens when there is no solution to the previous eq.
  {
    printf("There is not a solution to find the intial dt");
    dt_initial = 0;
  }
  return dt_initial;
}

#endif
