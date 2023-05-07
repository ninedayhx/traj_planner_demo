/**
 * @file map_server_node.cpp
 * @author ninedayhx (1170535490@qq.com)
 * @brief
 * @version 0.1
 * @date 2023-05-04
 *
 * @copyright Copyright (c) 2023
 *
 */
// c++
#include <iostream>
// ros
#include <ros/ros.h>
#include <ros/package.h>
#include <visualization_msgs/Marker.h>
// pcl
#include <pcl/io/pcd_io.h>
#include <pcl/point_types.h>
#include <pcl/conversions.h>
#include <pcl/point_cloud.h>
#include <pcl_conversions/pcl_conversions.h>
#include <sensor_msgs/PointCloud2.h>
// grid map
#include <grid_map_core/GridMap.hpp>
#include <grid_map_ros/GridMapRosConverter.hpp>
#include "../include/GridMapPclLoader.hpp"
#include "../include/helpers.hpp"

using namespace std;
namespace gm = ::grid_map::grid_map_pcl;

// 定义点云类型
typedef pcl::PointXYZ PointT;
typedef pcl::PointCloud<PointT> PointCloud;

int main(int argc, char **argv)
{
    // ros init
    ros::init(argc, argv, "map_node");
    ros::NodeHandle nh;
    ros::Publisher grid_map_pub = nh.advertise<grid_map_msgs::GridMap>("grid_map_demo", 1, true);
    ros::Publisher map_cloud_pub = nh.advertise<sensor_msgs::PointCloud2>("map_cloud", 1, true);

    // read pcd file
    PointCloud::Ptr map_cloud(new PointCloud);
    if (pcl::io::loadPCDFile(ros::package::getPath("traj_plan") + "/map/transform.pcd", *map_cloud) < 0)
    {
        ROS_ERROR("Couldn't read pcd file");
    }

    sensor_msgs::PointCloud2 cloud_msg;

    ros::Rate loop_rate(1);

    grid_map::GridMapPclLoader gridMapPclLoader;
    string pathToCloud = ros::package::getPath("traj_plan") + "/map/transform.pcd";
    gridMapPclLoader.loadParameters(gm::getParameterPath(nh));
    gridMapPclLoader.loadCloudFromPcdFile(pathToCloud);

    gm::processPointcloud(&gridMapPclLoader, nh);
    grid_map::GridMap gridMap = gridMapPclLoader.getGridMap();
    gridMap.setFrameId("map");

    grid_map_msgs::GridMap map_msg;
    grid_map::GridMapRosConverter::toMessage(gridMap, map_msg);

    grid_map_pub.publish(map_msg);

    // while (ros::ok())
    // {
    ROS_INFO("Published: map_msg");
    pcl::toROSMsg(*map_cloud, cloud_msg);
    cloud_msg.header.stamp = ros::Time::now();
    cloud_msg.header.frame_id = "map";

    map_cloud_pub.publish(cloud_msg);
    //     // 处理回调函数
    ros::spin();
    //     // 按照循环频率休眠
    //     loop_rate.sleep();
    // }
    return 0;
}