#include "ns3/core-module.h"
#include "ns3/network-module.h"
#include "ns3/point-to-point-module.h"
#include "ns3/internet-module.h"
#include "ns3/point-to-point-dumbbell.h"

using namespace ns3;

int main (int argc, char *argv[])
{
  // Set up logging

  // Set up the point-to-point links and devices
  PointToPointHelper p2p;
  p2p.SetDeviceAttribute ("DataRate", StringValue ("5Mbps"));
  p2p.SetChannelAttribute ("Delay", StringValue ("2ms"));

  // Dumbbell 1
  PointToPointDumbbellHelper dumbbell1(3, p2p, 3, p2p, p2p);
  
  // Dumbbell 2
  PointToPointDumbbellHelper dumbbell2(2, p2p, 2, p2p, p2p);
  
  // Dumbbell 3
  PointToPointDumbbellHelper dumbbell3(3, p2p, 3, p2p, p2p);

  // Connect the right side of dumbbell 1 to the left side of dumbbell 2
  p2p.Install(dumbbell1.GetRight(1), dumbbell2.GetLeft(0));

  // Connect the right side of dumbbell 2 to the left side of dumbbell 3
  p2p.Install(dumbbell2.GetRight(1), dumbbell3.GetLeft(0));

  // Set up the Internet stacks
  InternetStackHelper stack;
  dumbbell1.InstallStack(stack);
  dumbbell2.InstallStack(stack);
  dumbbell3.InstallStack(stack);

  // Set up IP addresses
  Ipv4AddressHelper leftIp = Ipv4AddressHelper("10.1.1.0", "255.255.255.0");
  Ipv4AddressHelper rightIp = Ipv4AddressHelper("10.2.1.0", "255.255.255.0");
  Ipv4AddressHelper routerIp = Ipv4AddressHelper("10.3.1.0", "255.255.255.0");

  dumbbell1.AssignIpv4Addresses(leftIp, rightIp, routerIp);
  leftIp.NewNetwork();
  rightIp.NewNetwork();
  routerIp.NewNetwork();

  dumbbell2.AssignIpv4Addresses(leftIp, rightIp, routerIp);
  leftIp.NewNetwork();
  rightIp.NewNetwork();
  routerIp.NewNetwork();

  dumbbell3.AssignIpv4Addresses(leftIp, rightIp, routerIp);

  // Populate global routing
  Ipv4GlobalRoutingHelper::PopulateRoutingTables();

  // Add your applications and simulation code here

  Simulator::Run ();
  Simulator::Destroy ();

  return 0;
}