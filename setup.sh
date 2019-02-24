#!/usr/bin/env bash

MYPWD=${PWD} 
cd ${MYPWD}/src/cvx/src/solvers/cvxgen_vel
make clean
make
cd ${MYPWD}/src/cvx/src/solvers/cvxgen_accel
make clean
make
cd ${MYPWD}/src/cvx/src/solvers/cvxgen_jerk
make clean
make
cd ${MYPWD}