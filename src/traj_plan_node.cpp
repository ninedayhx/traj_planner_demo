/**
 * @file trajrctory_visualization.cpp
 * @author ninedayhx (1170535490@qq.com)
 * @brief
 * @version 0.1
 * @date 2023-05-03
 *
 * @copyright Copyright (c) 2023
 *
 */
// c++
#include <iostream>
#include <vector>
// ros
#include <ros/ros.h>
#include <nav_msgs/Path.h>
#include <geometry_msgs/PointStamped.h>
// pcl
#include <pcl/io/pcd_io.h>
#include <pcl/point_types.h>
#include <pcl/conversions.h>
#include <pcl/point_cloud.h>
#include <pcl_conversions/pcl_conversions.h>
#include <sensor_msgs/PointCloud2.h>
// user: minimum_traj
#include "../include/minimum_traj.hpp"

using namespace std;

typedef pcl::PointXYZ PointT;
typedef pcl::PointCloud<PointT> PointCloud;

int main(int argc, char **argv)
{
    ros::init(argc, argv, "traj_plan");
    ros::NodeHandle nh;

    ros::Publisher path_pub = nh.advertise<nav_msgs::Path>("trajectory", 1000, true);
    ros::Publisher point_pub = nh.advertise<sensor_msgs::PointCloud2>("wawypoints", 1000, true);

    nav_msgs::Path path_msg;
    sensor_msgs::PointCloud2 waypoints;

    Eigen::MatrixXd postion3d(3, 5);
    Eigen::MatrixXd dstart3d(4, 3), dend3d(4, 3);
    Eigen::VectorXd Time(4);
    postion3d << 0.5, 0.3, 0.1, 0.2, 1.2,
        4.5, 1.5, -0.5, -2.3, -4.5,
        1.3, 1.2, 1.1, 1.1, 1.1;

    dstart3d << 0, 0, 0,
        0, 0, 0,
        0, 0, 0,
        0, 0, 0;
    dend3d << 0, 0, 0,
        0, 0, 0,
        0, 0, 0,
        0, 0, 0;
    Time << 1, 1, 1, 1;

    PointT point_tmp;
    PointCloud cloud_tmp;

    for (int k = 0; k < postion3d.cols(); k++)
    {
        point_tmp.x = postion3d(0, k);
        point_tmp.y = postion3d(1, k);
        point_tmp.z = postion3d(2, k);
        cloud_tmp.push_back(point_tmp);
    }

    pcl::toROSMsg(cloud_tmp, waypoints);
    waypoints.header.frame_id = "map";
    waypoints.header.stamp = ros::Time::now();
    point_pub.publish(waypoints);

    minimum_traj jerk(JERK, 3, 5, 3, postion3d, dstart3d, dend3d, Time);
    Eigen::MatrixXd traj = jerk.Cal_minimum_traj(Time, 3);

    path_msg.header.frame_id = "map";

    path_msg.header.stamp = ros::Time::now();

    for (int i = 0; i < Time.size(); i++)
    {
        for (int j = 0; j < 100; j++)
        {
            geometry_msgs::PoseStamped pose;
            pose.pose.position.x = traj(0, i * 100 + j);
            pose.pose.position.y = traj(1, i * 100 + j);
            pose.pose.position.z = traj(2, i * 100 + j);
            pose.pose.orientation.w = 1.0;
            path_msg.poses.push_back(pose);
        }
    }

    path_pub.publish(path_msg);
    ros::spin();

    return 0;
}
