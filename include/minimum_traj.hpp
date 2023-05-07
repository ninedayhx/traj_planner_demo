/**
 * @file minimum_traj.hpp
 * @author ninedayhx (1170535490@qq.com)
 * @brief
 * @version 0.1
 * @date 2023-04-15
 *
 * @copyright Copyright (c) 2023
 *
 */
#ifndef __MINI_TRAJ_H
#define __MINI_TRAJ_H
#include <iostream>
#include <fstream>
#include "eigen3/Eigen/Eigen"
#include <cmath>
#include <vector>

using namespace std;

#define JERK 3
#define SNAP 4

class minimum_traj
{
private:
    /* data */

    // store the points data provided by people,
    //
    //      the coordinate of start point,waypoints and end points
    //      the velocity of start points and end point
    //      the acceleration of start point and end point
    //
    // eg: a trajectory with n segments i.e. n+1 points; n-1 waypoints
    //    the structure of D_fixed will be a (2*3+n)*3 matrix
    //    x    y    z
    //    x0   y0   z0
    //    vx0  vy0  vz0
    //    ax0  ay0  az0
    //    xn   yn   zn
    //    vxn  vyn  vzn
    //    axn  ayn  azn
    //    x1   y1   z1
    //    x2   y2   z2
    //        .....
    //    xn-1 yn-1 zn-1
    //    which contains all the data people defined
    Eigen::MatrixXd D_fixed;
    // close form solution optimal result
    Eigen::MatrixXd D_P_optimal;
    // the solved waypoints of all trajectory
    Eigen::MatrixXd D_total;
    Eigen::MatrixXd D_total_selected;
    // the mapping matrix of all trajectory
    Eigen::MatrixXd A_total;
    Eigen::MatrixXd A_one;
    Eigen::MatrixXd A_one_t;
    // the hessis matrix of all trajectory
    Eigen::MatrixXd Q_total;
    Eigen::MatrixXd Q_one;
    Eigen::MatrixXd Q_one_t;
    // the Transport of select matrix to devide the waypoints data into fixed and free
    Eigen::MatrixXd C_select_T;
    Eigen::MatrixXd R;
    Eigen::MatrixXd R_FF;
    Eigen::MatrixXd R_FP;
    Eigen::MatrixXd R_PF;
    Eigen::MatrixXd R_PP;

    // unsigned int poly_order;
    unsigned int poly_coff_num;
    unsigned int all_poly_coff_num;
    // unsigned int physical_num;
    unsigned int points_num;
    unsigned int seg_num;
    unsigned int fixed_coff_num;
    unsigned int all_coff_num;

    /* function */

    int Fac(int x);

public:
    /* data */
    Eigen::MatrixXd Poly_coff_total;
    /* function */
    minimum_traj(
        unsigned int minimum_type,
        unsigned int dim,
        unsigned int poly_order,
        unsigned int physical_num,
        const Eigen::MatrixXd &Pos,
        const Eigen::MatrixXd &dStart,
        const Eigen::MatrixXd &dEnd,
        const Eigen::VectorXd &Time);

    ~minimum_traj();

    bool Cal_C_select_T(Eigen::MatrixXd &C_T, unsigned int physical_num);

    Eigen::MatrixXd Cal_minimum_traj(Eigen::VectorXd &Time, int dim);

    Eigen::MatrixXd Get_Poly_coff_total();
};

#endif