#include "../include/minimum_traj.hpp"

/**
 * @brief Construct a new minimum traj::minimum traj object
 *
 * @param minimum_type
 * @param dim
 * @param poly_order
 * @param physical_num
 * @param Pos
 * @param dStart
 * @param dEnd
 * @param Time
 */
minimum_traj::minimum_traj(
    unsigned int minimum_type,
    unsigned int dim,
    unsigned int poly_order,
    unsigned int physical_num,
    const Eigen::MatrixXd &Pos,
    const Eigen::MatrixXd &dStart,
    const Eigen::MatrixXd &dEnd,
    const Eigen::VectorXd &Time)
{
    // pos is a 3*n matrix D_fixed is a n*3 matrix
    // for convenience， transform the pos to n*3
    Eigen::MatrixXd Pos_tmp = Pos.transpose();

    poly_coff_num = poly_order + 1;

    points_num = Time.size() + 1;
    seg_num = Time.size();

    fixed_coff_num = physical_num * 2 + (points_num - 2) * 2;
    all_coff_num = physical_num * (2 + 2 * (points_num - 2));

    D_fixed = Eigen::MatrixXd::Zero(fixed_coff_num, dim);
    D_P_optimal = Eigen::MatrixXd::Zero(all_coff_num - fixed_coff_num, dim);
    D_total = Eigen::MatrixXd::Zero(all_coff_num, dim);
    D_total_selected = Eigen::MatrixXd::Zero(all_coff_num, dim);

    A_total = Eigen::MatrixXd::Zero(seg_num * physical_num * 2, seg_num * poly_coff_num);
    A_one = Eigen::MatrixXd::Zero(physical_num * 2, poly_coff_num);
    A_one_t = Eigen::MatrixXd::Zero(physical_num * 2, poly_coff_num);

    Q_total = Eigen::MatrixXd::Zero(seg_num * poly_coff_num, seg_num * poly_coff_num);
    Q_one = Eigen::MatrixXd::Zero(poly_coff_num, poly_coff_num);
    Q_one_t = Eigen::MatrixXd::Zero(poly_coff_num, poly_coff_num);

    C_select_T = Eigen::MatrixXd::Zero(all_coff_num, D_total.rows());

    R = Eigen::MatrixXd::Zero(all_coff_num, all_coff_num);

    // dStart(physical_num-1,dim)
    // dEnd(physical_num-1,dim)
    // the solve of jerk and snap is meanless
    // if (physical_num > 3) the jerk and snap will be zero
    D_total.row(0) = Pos_tmp.row(0);
    D_total.row(all_coff_num - physical_num) = Pos_tmp.row(points_num - 1);
    for (int i = 0; i < points_num - 2; i++)
    {
        D_total.row(physical_num + i * physical_num * 2) = Pos_tmp.row(i + 1);
        D_total.row(physical_num + i * physical_num * 2 + physical_num) = Pos_tmp.row(i + 1);
    }
    for (int i = 0; i < physical_num - 1; i++)
    {
        D_total.row(1 + i) = dStart.row(i);
        D_total.row(all_coff_num - physical_num + 1 + i) = dEnd.row(i);
    }

    // init A_total
    // A_total.resize(seg_num * 6, seg_num * 6);
    // caluate A by time segments
    for (int row = 0; row < physical_num; row++)
    {
        for (int col = 0; col < poly_coff_num; col++)
        {
            if (row == col)
            {
                A_one(row, col) = Fac(row);
            }
            if (row <= col)
            {
                A_one(row + physical_num, col) = Fac(col) / Fac(col - row);
            }
        }
    }

    for (int M = 0; M < seg_num; M++)
    {
        for (int i = 0; i < physical_num; i++)
        {
            for (int j = 0; j < poly_coff_num; j++)
            {
                int order = j - i;
                if (order < 0)
                {
                    order = 0;
                }
                A_one_t(i, j) = A_one(i, j) * pow(0, order);
                A_one_t(i + physical_num, j) = A_one(i + physical_num, j) * pow(Time(M), order);
            }
        }
        A_total.block(M * poly_coff_num, M * poly_coff_num, poly_coff_num, poly_coff_num) = A_one_t;
        A_one_t = Eigen::MatrixXd::Zero(poly_coff_num, poly_coff_num);
    }

    // init Q
    unsigned int k = minimum_type;

    for (int m = k; m < poly_coff_num; m++)
    {
        for (int n = k; n < poly_coff_num; n++)
        {
            Q_one(m, n) = Fac(m) * Fac(n) / (Fac(m - k) * Fac(n - k) * (m + n - 2 * k + 1));
        }
    }

    for (int M = 0; M < seg_num; M++)
    {
        for (int m = k; m < poly_coff_num; m++)
        {
            for (int n = k; n < poly_coff_num; n++)
            {
                Q_one_t(m, n) = Q_one(m, n) * pow(Time(M), m + n - 2 * k + 1);
            }
        }
        Q_total.block(M * poly_coff_num, M * poly_coff_num, poly_coff_num, poly_coff_num) = Q_one_t;
        Q_one_t = Eigen::MatrixXd::Zero(poly_coff_num, poly_coff_num);
    }

    // init select matrix
    if (!Cal_C_select_T(C_select_T, physical_num))
    {
        return;
    }

    D_total_selected = C_select_T * D_total;
    D_fixed = D_total_selected.block(0, 0, fixed_coff_num, dim);

    R = C_select_T.transpose() * A_total.transpose().inverse() * Q_total * A_total.inverse() * C_select_T;

    R_FF = R.block(0, 0, fixed_coff_num, fixed_coff_num);
    R_FP = R.block(0, fixed_coff_num, fixed_coff_num, R.cols() - fixed_coff_num);
    R_PF = R.block(fixed_coff_num, 0, R.rows() - fixed_coff_num, fixed_coff_num);
    R_PP = R.block(fixed_coff_num, fixed_coff_num, R.rows() - fixed_coff_num, R.cols() - fixed_coff_num);

    D_P_optimal = -R_PP.inverse() * R_FP.transpose() * D_fixed;

    cout << "R size:" << R.rows() << "," << R.cols() << endl;
    cout << C_select_T.transpose().rows() << "*" << C_select_T.transpose().cols() << "|"
         << A_total.transpose().inverse().rows() << "*" << A_total.transpose().inverse().cols() << "|"
         << Q_total.rows() << "*" << Q_total.cols() << "|"
         << A_total.inverse().rows() << "*" << A_total.inverse().cols() << "|"
         << C_select_T.rows() << "*" << C_select_T.cols() << endl;

    // D_total_selected.block(0, 0, fixed_coff_num, 3) = D_fixed;
    D_total_selected.block(fixed_coff_num, 0, all_coff_num - fixed_coff_num, dim) = D_P_optimal;

    Poly_coff_total = A_total.inverse() * C_select_T.transpose() * D_total_selected;
    // std::ofstream fout("poly_coff_all.csv", std::ios::binary);
    // fout << Poly_coff_total << std::endl;
    // fout.flush();
    cout << "Poly_coff_total" << endl
         << Poly_coff_total << endl;
}

minimum_traj::~minimum_traj()
{
    cout << "delete" << endl;
}

bool minimum_traj::Cal_C_select_T(Eigen::MatrixXd &C_T, unsigned int phy_num)
{
    int tmp[10] = {};
    cout << "test" << endl;
    C_T = Eigen::MatrixXd::Zero(D_total.rows(), D_total.rows());
    for (int i = 0; i < all_coff_num; i++)
    {

        // start point vel acc
        if (i < phy_num)
        {
            C_T(i, i) = 1;
            continue;
        }
        // end point vel acc
        if ((i >= phy_num) && (i < 2 * phy_num))
        {
            C_T(i, (C_T.cols() - phy_num) + (i - phy_num)) = 1;
            continue;
        }
        // waypoints
        if ((i >= 2 * phy_num) && (i < fixed_coff_num))
        {
            tmp[0]++;
            C_T(i, tmp[0] * phy_num) = 1;
            continue;
        }

        // waypoint derivate
        if ((i >= fixed_coff_num) && (i < all_coff_num))
        {
            // 用来看是第几阶微分 1-phy_num
            int der_cnt = (i - fixed_coff_num) / (2 * (points_num - 2)) + 1;
            // 用来看是第几个点 0-航点数-1
            int pt_cnt = (i - fixed_coff_num) % (2 * (points_num - 2));
            // tmp[der_cnt]++;
            //  起始偏置+第几个点+选择微分
            C_T(i, phy_num + pt_cnt * phy_num + der_cnt) = 1;
            continue;
        }
    }
    cout << C_T << endl;
    if (C_T.determinant() == 0)
    {
        cout << "C_select is error";
        return false;
    }
    if (C_T.inverse() != C_T.transpose())
    {
        cout << "C_select is error";
        return false;
    }
    return true;
}

int minimum_traj::Fac(int x)
{
    int fac = 1;
    for (int i = x; i > 0; i--)
        fac = fac * i;
    return fac;
}

Eigen::MatrixXd minimum_traj::Cal_minimum_traj(Eigen::VectorXd &Time, int dim)
{
    unsigned int seg_num = Time.size();
    Eigen::MatrixXd traj_res = Eigen::MatrixXd::Zero(dim, 100 * seg_num);
    for (int i = 0; i < seg_num; i++)
    {
        for (int j = 0; j < 100; j++)
        {
            double t = j * 0.01;
            for (int k = 0; k < 6; k++)
            {
                for (int d = 0; d < dim; d++)
                {
                    traj_res(d, i * 100 + j) += Poly_coff_total(i * 6 + k, d) * pow(t, k);
                }
            }
        }
    }
    // std::ofstream fout("matrixTest.csv", std::ios::binary);
    // fout << traj_res << std::endl;
    // fout.flush();
    return traj_res;
}

Eigen::MatrixXd minimum_traj::Get_Poly_coff_total()
{
    return Poly_coff_total;
}