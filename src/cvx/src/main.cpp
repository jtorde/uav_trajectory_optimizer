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

#include "solvers/solvers.hpp"

int main(int argc, char **argv)
{
  // Solvers objects
  Solver<VEL> solver_vel;
  Solver<ACCEL> solver_accel;
  Solver<JERK> solver_jerk;

  double v_max = 5;
  double a_max = 10;
  double j_max = 20;
  double max_values[3] = { v_max, a_max, j_max };

  // Example 1: Example when input=Velocity
  printf("\n\n****** Example 1: Input=Velocity *****\n");
  double x0_ex1[6] = { 1, 3, 4 };  // Initial State: x,y,z
  double xf_ex1[6] = { 2, 8, 3 };  // Final State: x,y,z
  solver_vel.setq(200000);         // Weight of the final cost
  solver_vel.set_max(max_values);

  solver_vel.set_x0(x0_ex1);
  solver_vel.set_xf(xf_ex1);
  solver_vel.genNewTraj();

  double **x_ex1;
  double **u_ex1;
  x_ex1 = solver_vel.getState();
  u_ex1 = solver_vel.getInput();

  printf("\nPositions:          Velocities: \n");
  printf("%0.2f %0.2f %0.2f      %0.2f %0.2f %0.2f \n", x0_ex1[0], x0_ex1[1], x0_ex1[2], u_ex1[0][0], u_ex1[0][1],
         u_ex1[0][2]);
  for (int i = 1; i < N_VEL; i++)
  {
    printf("%0.2f %0.2f %0.2f      %0.2f %0.2f %0.2f \n", x_ex1[i][0], x_ex1[i][1], x_ex1[i][2], u_ex1[i][0],
           u_ex1[i][1], u_ex1[i][2]);
  }

  printf("Cost= %0.2f\n", solver_vel.getCost());

  // Example 2: Example when input=Acceleration
  printf("\n\n****** Example 2: Input=Acceleration *****\n");
  double x0_ex2[6] = { 1, 3, 4, 5, 7, 3 };  // Initial State: x,y,z,vx,vy,vz
  double xf_ex2[6] = { 2, 8, 3, 4, 9, 2 };  // Final State: x,y,z,vx,vy,vz
  solver_accel.setq(200000);                // Weight of the final cost
  solver_accel.set_max(max_values);

  solver_accel.set_x0(x0_ex2);
  solver_accel.set_xf(xf_ex2);
  solver_accel.genNewTraj();

  double **x_ex2;
  double **u_ex2;
  x_ex2 = solver_accel.getState();
  u_ex2 = solver_accel.getInput();

  printf("\nPositions:          Velocities:          Accelerations: \n");
  printf("%0.2f %0.2f %0.2f      %0.2f %0.2f %0.2f      %0.2f %0.2f %0.2f \n", x0_ex2[0], x0_ex2[1], x0_ex2[2],
         x0_ex2[3], x0_ex2[4], x0_ex2[5], u_ex2[0][0], u_ex2[0][1], u_ex2[0][2]);
  for (int i = 1; i < N_ACCEL; i++)
  {
    printf("%0.2f %0.2f %0.2f      %0.2f %0.2f %0.2f      %0.2f %0.2f %0.2f\n", x_ex2[i][0], x_ex2[i][1], x_ex2[i][2],
           x_ex2[i][3], x_ex2[i][4], x_ex2[i][5], u_ex2[i][0], u_ex2[i][1], u_ex2[i][2]);
  }

  printf("\nCost= %0.2f\n", solver_accel.getCost());

  // Example 3: Example when input=Jerk
  printf("\n\n****** Example 3: Input=Jerk *****\n");
  double x0_ex3[9] = { 1, 3, 4, 5, 7, 3, 6, 7, 3 };  // Initial State: x,y,z.vx,vy,vz,ax,ay,az
  double xf_ex3[9] = { 2, 8, 3, 4, 9, 2, 2, 4, 1 };  // Final State: x,y,z.vx,vy,vz,ax,ay,az
  solver_jerk.setq(200000);                          // Weight of the final cost
  solver_jerk.set_max(max_values);

  solver_jerk.set_x0(x0_ex3);
  solver_jerk.set_xf(xf_ex3);
  solver_jerk.genNewTraj();

  double **x_ex3;
  double **u_ex3;
  x_ex3 = solver_jerk.getState();
  u_ex3 = solver_jerk.getInput();

  printf("\nPositions:          Velocities:          Accelerations:    Jerks:\n");
  printf("%0.2f %0.2f %0.2f      %0.2f %0.2f %0.2f      %0.2f %0.2f %0.2f     %0.2f %0.2f %0.2f\n", x0_ex3[0],
         x0_ex3[1], x0_ex3[2], x0_ex3[3], x0_ex3[4], x0_ex3[5], x0_ex3[6], x0_ex3[7], x0_ex3[8], u_ex3[0][0],
         u_ex3[0][1], u_ex3[0][2]);
  for (int i = 1; i < N_JERK; i++)
  {
    printf("%0.2f %0.2f %0.2f      %0.2f %0.2f %0.2f      %0.2f %0.2f %0.2f     %0.2f %0.2f %0.2f\n", x_ex3[i][0],
           x_ex3[i][1], x_ex3[i][2], x_ex3[i][3], x_ex3[i][4], x_ex3[i][5], x_ex3[i][6], x_ex3[i][7], x_ex3[i][8],
           u_ex3[i][0], u_ex3[i][1], u_ex3[i][2]);
  }
  printf("\nCost= %0.2f\n", solver_jerk.getCost());

  return 0;
}
