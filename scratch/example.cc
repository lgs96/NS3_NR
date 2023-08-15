#include "ns3/core-module.h"
#include "ns3/network-module.h"
#include "ns3/internet-module.h"
#include "ns3/point-to-point-module.h"
#include "ns3/point-to-point-dumbbell.h"
#include "ns3/applications-module.h"
#include "ns3/ipv4-global-routing-helper.h"
#include "ns3/log.h"

#include <fstream>
#include <iostream>

using namespace ns3;

NS_LOG_COMPONENT_DEFINE("OpcSimulation");

class ObjectHeader : public Header 
{
public:
  ObjectHeader () : m_index (0), m_size (0) {}

  static TypeId GetTypeId (void) {
    static TypeId tid = TypeId ("ObjectHeader")
      .SetParent<Header> ()
      .AddConstructor<ObjectHeader> ();
    return tid;
  }

  virtual TypeId GetInstanceTypeId (void) const { return GetTypeId (); }
  virtual void Print (std::ostream &os) const { os << "Index: " << m_index << ", Size: " << m_size; }
  virtual uint32_t GetSerializedSize (void) const { return 8; }
  virtual void Serialize (Buffer::Iterator start) const {
    start.WriteU32 (m_index);
    start.WriteU32 (m_size);
  }
  virtual uint32_t Deserialize (Buffer::Iterator start) {
    m_index = start.ReadU32 ();
    m_size = start.ReadU32 ();
    return 8;
  }

  void SetIndex (uint32_t index) { m_index = index; }
  uint32_t GetIndex () const { return m_index; }

  void SetSize (uint32_t size) { m_size = size; }
  uint32_t GetSize () const { return m_size; }

private:
  uint32_t m_index;
  uint32_t m_size;
};

class Logger {
public:
  static Logger& Instance(const std::string &filename = "") {
    static Logger instance;
    if (!filename.empty()) {
      instance.filename_ = filename;
    }
    return instance;
  }

  void WriteCSVHeader() {
    std::ofstream csvFile;
    csvFile.open(filename_);
    csvFile << "Number\tTransmission time\tComputing time\tTotal time\tTraversed node num\n";
    csvFile.close();
  }

  void LogCompletedObject(double txTime, double computeTime, double totalTime, int traversedNodeNum) {
    std::ofstream csvFile;
    csvFile.open(filename_, std::ios_base::app); // Append to the file
    number += 1;
    csvFile << number << "\t" << txTime << "\t" << computeTime << "\t" << totalTime << "\t" << traversedNodeNum << "\n";
    csvFile.close();
  }

private:
  std::string filename_;
  uint32_t number = 0;

  // Private constructor to ensure only one instance is created
  Logger() {}
  Logger(const Logger&) = delete;
  Logger& operator=(const Logger&) = delete;
};

class ObjectTag : public Tag
{
public:
  enum Log { GENERATED, TRANSMITTED, COMPUTED, COMPLETED };

  ObjectTag () : m_eventTime(0), m_nodesTraversed(0), m_transmissionTime(0), m_computingTime(0) {}

  static TypeId GetTypeId (void) {
    static TypeId tid = TypeId ("ObjectTag")
      .SetParent<Tag> ()
      .AddConstructor<ObjectTag> ();
    return tid;
  }

  virtual TypeId GetInstanceTypeId (void) const {
    return GetTypeId ();
  }

  virtual uint32_t GetSerializedSize (void) const {
    return 4 * sizeof(double); // Four double variables
  }

  virtual void Serialize (TagBuffer i) const {
    i.WriteDouble(m_eventTime);
    i.WriteDouble(m_nodesTraversed);
    i.WriteDouble(m_transmissionTime);
    i.WriteDouble(m_computingTime);
  }

  virtual void Deserialize (TagBuffer i) {
    m_eventTime = i.ReadDouble();
    m_nodesTraversed = i.ReadDouble();
    m_transmissionTime = i.ReadDouble();
    m_computingTime = i.ReadDouble();
  }

  virtual void Print (std::ostream &os) const {
    os << "Event Time: " << m_eventTime << ", Nodes Traversed: " << m_nodesTraversed << ", Transmission Time: " << m_transmissionTime << ", Computing Time: " << m_computingTime;
  }

  void SetEventTime (double time) { m_eventTime = time; }
  double GetEventTime () const { return m_eventTime; }

  void SetNodesTraversed (double nodes) { m_nodesTraversed = nodes; }
  double GetNodesTraversed () const { return m_nodesTraversed; }

  void SetTransmissionTime (double time) { m_transmissionTime = time; }
  double GetTransmissionTime () const { return m_transmissionTime; }

  void SetComputingTime (double time) { m_computingTime = time; }
  double GetComputingTime () const { return m_computingTime; }

  void LogEvent(Log eventType) {
    double currentTime = Simulator::Now().GetSeconds();
    switch (eventType) {
      case GENERATED:
        NS_LOG_DEBUG("Event GENERATED at time: " << m_eventTime << "s");
        break;
      case TRANSMITTED:
        m_transmissionTime = currentTime - m_eventTime;
        NS_LOG_DEBUG("Event TRANSMITTED at time: " << currentTime << "s, Transmission Time: " << m_transmissionTime << "s");
        break;
      case COMPUTED:
        m_computingTime = currentTime - m_eventTime;
        m_nodesTraversed += 1;
        NS_LOG_DEBUG("Event COMPUTED at time: " << currentTime << "s, Computing Time: " << m_computingTime << "s");
        break;
      case COMPLETED:
        m_computingTime = currentTime - m_eventTime;
        m_nodesTraversed += 1;
        Logger::Instance().LogCompletedObject(m_transmissionTime, m_computingTime, m_transmissionTime + m_computingTime, m_nodesTraversed);
        NS_LOG_DEBUG("Event COMPLETED at time: " << currentTime << "s");
        break;
    }
    m_eventTime = currentTime;
  }

private:
  double m_eventTime;
  double m_nodesTraversed;
  double m_transmissionTime;
  double m_computingTime;
};

class ReceiveSendApplication : public Application 
{
public:
  enum Mode { RECEIVE_ONLY, TRANSMIT_ONLY, RECEIVE_AND_TRANSMIT };

  ReceiveSendApplication (Mode mode, uint32_t objectSize, Time transmitInterval, double exitProb, double compressRatio) : m_mode(mode), m_transmitInterval(transmitInterval), 
                m_currentObjectSize(objectSize), m_currentObjectReceivedSize(0), m_maxPacketSize(512), m_earlyExitProb(exitProb), m_compressRatio(compressRatio) {}
  virtual ~ReceiveSendApplication() {}

  void Setup (Ptr<Socket> recvSocket, Ptr<Socket> sendSocket, Address sendAddress, uint32_t packetSize, Time interval) {
    m_recvSocket = recvSocket;
    m_sendSocket = sendSocket;
    m_sendAddress = sendAddress;
    m_packetSize = packetSize;
    m_interval = interval;
  }

private:
   virtual void StartApplication (void) {
    m_running = true;

    if (m_mode != TRANSMIT_ONLY) {
        m_recvSocket->Bind (InetSocketAddress (Ipv4Address::GetAny (), 8080));
        m_recvSocket->Listen ();
        m_recvSocket->SetAcceptCallback(
            MakeNullCallback<bool, Ptr<Socket>, const Address&>(),
            MakeCallback (&ReceiveSendApplication::HandleAccept, this)
        );
    }

    if (m_mode != RECEIVE_ONLY) {
        InetSocketAddress sendSocketAddress = InetSocketAddress::ConvertFrom(m_sendAddress);
        Ipv4Address sendIpv4Address = sendSocketAddress.GetIpv4();
        NS_LOG_UNCOND(Simulator::Now().GetSeconds()<<" Connection request to "<< sendIpv4Address <<" from "<< m_sendSocket->GetNode()->GetObject<Ipv4>()->GetAddress(1, 0).GetLocal());
        m_sendSocket->Connect(m_sendAddress);
        m_sendSocket->SetConnectCallback (MakeCallback (&ReceiveSendApplication::ConnectionSucceededCallback,
                                                    this),
                                        MakeCallback (&ReceiveSendApplication::ConnectionFailedCallback,
                                                    this));
    }
    // rest of the method

    if (m_mode == TRANSMIT_ONLY) {
      m_transmitEvent = Simulator::ScheduleNow (&ReceiveSendApplication::TransmitPeriodically, this);
    }
  }

  bool HandleConnectionRequest (Ptr<Socket> s, const Address& from) {
    InetSocketAddress remoteSocketAddress = InetSocketAddress::ConvertFrom(from);
    Ipv4Address remoteIpv4Address = remoteSocketAddress.GetIpv4();
    
    m_remoteAddress = remoteIpv4Address;

    NS_LOG_UNCOND (Simulator::Now().GetSeconds()<<" Received a connection request from " << remoteIpv4Address<< " to "<< s->GetNode()->GetObject<Ipv4>()->GetAddress(1, 0).GetLocal());

    return true; // Accept the connection request
  }


  void HandleAccept (Ptr<Socket> s, const Address& from) {
    s->SetRecvCallback(MakeCallback(&ReceiveSendApplication::ReceivePacket, this));
    m_connectedSockets.push_back(s);

    NS_LOG_UNCOND ("Accepted connection from: " << InetSocketAddress::ConvertFrom (from).GetIpv4 ());
  }

  void ConnectionSucceededCallback (Ptr<Socket> socket)
    {
        NS_LOG_UNCOND (Simulator::Now().GetSeconds() <<" " << socket->GetNode()->GetObject<Ipv4>()->GetAddress(1, 0).GetLocal() << " Connection done");
    }

   void ConnectionFailedCallback (Ptr<Socket> socket)
   {
    NS_LOG_UNCOND (Simulator::Now().GetSeconds() <<" " << socket->GetNode()->GetObject<Ipv4>()->GetAddress(1, 0).GetLocal() << " Connection failed");
   }


  virtual void StopApplication (void) {
    m_running = false;
    if (m_sendEvent.IsRunning ())
      Simulator::Cancel (m_sendEvent);

    if (m_recvSocket)
      m_recvSocket->Close ();

    for (Ptr<Socket> socket : m_connectedSockets) {
      socket->Close ();
    }

    if (m_sendSocket)
      m_sendSocket->Close ();
  }

  void TransmitPeriodically () {
    if (m_running) {
      if (m_mode == TRANSMIT_ONLY) {
        SendObject (m_currentObjectSize); // You can set the size of objects to send
      }

      m_transmitEvent = Simulator::Schedule (m_transmitInterval, &ReceiveSendApplication::TransmitPeriodically, this);
    }
  }

  void ReceivePacket (Ptr<Socket> socket) {
    Ptr<Packet> packet;
    ObjectHeader header;
    while ((packet = socket->Recv ()))
    {
      packet->RemoveHeader (header);
      
      // Check if this is the start of a new object
      if (header.GetIndex() == 0) {
        m_currentObjectSize = header.GetSize();
        m_currentObjectReceivedSize = 0;
      }

      m_currentObjectReceivedSize += packet->GetSize();
      //NS_LOG_UNCOND (Simulator::Now().GetSeconds()<<" Received entire object at node " << socket->GetNode()->GetId()<<" "<< m_remoteAddress<< " "<< m_sendAddress);
      // Check if we have received the entire object
      if (m_currentObjectReceivedSize == m_currentObjectSize) {
        Ipv4Address fromNodeIp = m_remoteAddress;
        uint32_t toNodeId = socket->GetNode()->GetId();
        Ipv4Address toNodeIp = socket->GetNode()->GetObject<Ipv4>()->GetAddress(1, 0).GetLocal(); // Assuming the IP is at interface 1
        NS_LOG_UNCOND (Simulator::Now().GetSeconds() << " Received entire object at node " << toNodeId << " (" << toNodeIp << ") from node " << fromNodeIp << " with size " << m_currentObjectSize);
   
        double randomValue = static_cast<double>(rand()) / RAND_MAX;

        if (m_mode == RECEIVE_AND_TRANSMIT) {
            // Transmitting an object of the same size to the next hierarchy
            NS_LOG_UNCOND (Simulator::Now().GetSeconds()<<" Relay received object to the next node");
            SendObject (m_currentObjectSize);
        }
        else {
            NS_LOG_UNCOND (Simulator::Now().GetSeconds()<<" Object inference is completed in the server node");
            return;
        }
      }
    }
  }

  void SendObject (uint32_t objectSize) {

    objectSize = objectSize*m_compressRatio;

    if (m_running) {
      uint32_t fromNodeId = m_sendSocket->GetNode()->GetId();
      Ipv4Address fromNodeIp = m_sendSocket->GetNode()->GetObject<Ipv4>()->GetAddress(1, 0).GetLocal(); // Assuming the IP is at interface 1
      Ipv4Address toNodeIp = InetSocketAddress::ConvertFrom(m_sendAddress).GetIpv4();
      uint32_t remainingSize = objectSize;
      uint32_t index = 0;

      double randomValue = static_cast<double>(rand()) / RAND_MAX;
      if (randomValue <= m_earlyExitProb) {
        NS_LOG_UNCOND ("Early exit for node " << fromNodeId);
        return;
      }

      // Divide the object into packets
      while (remainingSize > 0) {
        uint32_t packetSize = std::min(remainingSize, m_maxPacketSize);
        Ptr<Packet> packet = Create<Packet> (packetSize);
        ObjectHeader header;
        header.SetIndex (index++);
        header.SetSize (objectSize);
        packet->AddHeader (header);
        m_sendSocket->Send (packet);

        remainingSize -= packetSize;
      }

      NS_LOG_UNCOND (Simulator::Now().GetSeconds() << " Sent entire object from node " << fromNodeId << " (" << fromNodeIp << ") to " << toNodeIp << " with size " << objectSize);
   }
  }


  Ptr<Socket> m_recvSocket;
  Ptr<Socket> m_sendSocket;
  Address m_sendAddress;
  uint32_t m_packetSize;
  Time m_interval;
  EventId m_sendEvent;
  bool m_running = false;

  uint32_t m_currentObjectSize;
  uint32_t m_currentObjectIndex;
  uint32_t m_currentObjectReceivedSize;
  uint32_t m_maxPacketSize;

  double m_earlyExitProb;
  double m_compressRatio;

  uint32_t m_parallelComputingCapacitiy;  // weak-scale parameter
  double m_sequentialComputingAbility; // strong-scale parameter --> unit: 

  Mode m_mode;
  EventId m_transmitEvent;
  Time m_transmitInterval;

  private:
    Ipv4Address m_remoteAddress;
    std::vector<Ptr<Socket>> m_connectedSockets;
};

int main (int argc, char *argv[])
{

  //LogComponentEnable("Ipv4GlobalRouting", LOG_LEVEL_ALL);
  //LogComponentEnable("Ipv4StaticRouting", LOG_LEVEL_ALL);
  //LogComponentEnable("Ipv4L3Protocol", LOG_LEVEL_ALL);

  // Parameters
  uint32_t numberOfServers = 1;
  uint32_t numberOfCentralNodes = 3;
  uint32_t numberOfLowerNodes = 2; // Per central node, therefore total 4 nodes in this example
  uint32_t packetSize = 100; // xx bytes
  Time interval = MilliSeconds(50); // yy milliseconds

  // Create Nodes
  NodeContainer serverNodes;
  serverNodes.Create (numberOfServers);

  NodeContainer centralNodes;

  NodeContainer lowerNodes;

   // Create PointToPoint links
  PointToPointHelper p2p;
  p2p.SetDeviceAttribute ("DataRate", StringValue ("5Mbps"));
  p2p.SetChannelAttribute ("Delay", StringValue ("2ms"));

  // First level dumbbells
  centralNodes.Create(numberOfCentralNodes);
  std::vector<PointToPointDumbbellHelper> firstLevelDumbbells;

  Ipv4AddressHelper leftIp, rightIp, routerIp;
  leftIp.SetBase("10.1.0.0", "255.255.255.0");
  rightIp.SetBase("10.2.0.0", "255.255.255.0");
  routerIp.SetBase("10.3.0.0", "255.255.255.0");
  
  // Install Internet Stack
  InternetStackHelper stack;
  stack.Install (serverNodes);
  //internet.Install (centralNodes);
  //internet.Install (lowerNodes);

  
Ipv4AddressHelper ipv4h;
Ipv4StaticRoutingHelper ipv4RoutingHelper;
ipv4h.SetBase("10.4.0.0", "255.255.255.0");

  
  for (uint32_t i = 0; i < numberOfCentralNodes; ++i) {

    PointToPointDumbbellHelper d(numberOfLowerNodes, p2p, 1, p2p, p2p);

    NS_LOG_UNCOND("Set up first-level dumbbell "<< i);
    

    for (uint32_t i = 0; i < d.LeftCount(); ++i)
    {
        stack.Install(d.GetLeft(i));
    }
    for (uint32_t i = 0; i < d.RightCount(); ++i)
    {
        stack.Install(d.GetRight(i));
    }
    stack.Install(d.GetLeft());
    stack.Install(d.GetRight());

    firstLevelDumbbells.push_back(d);
    centralNodes.Add(d.GetRight(0));
    d.AssignIpv4Addresses(leftIp, rightIp, routerIp);


    for (uint32_t i = 0; i < d.RightCount(); ++i)
    {
        Ptr<Node> node = d.GetRight(i); // Replace with the appropriate node container and index
        Ptr<Ipv4> ipv4 = node->GetObject<Ipv4>();
        Ipv4Address addr = ipv4->GetAddress(1,0).GetLocal(); // Get the IP address of the first interface
        std::cout << "IP address of node is " << addr << std::endl;
    }

    
    std::string leftBase = "10.1." + std::to_string((i+1)*10) + ".0";
    std::string rightBase = "10.2." + std::to_string((i+1)*10) + ".0";
    std::string routerBase = "10.3." + std::to_string((i+1)*10) + ".0";

    leftIp.SetBase(leftBase.c_str(), "255.255.255.0");
    rightIp.SetBase(rightBase.c_str(), "255.255.255.0");
    routerIp.SetBase(routerBase.c_str(), "255.255.255.0");

    NS_LOG_UNCOND(i << "[First] Size of firstlevelDumbell: "<< firstLevelDumbbells[i].LeftCount()<< " "<< firstLevelDumbbells[i].RightCount());

    // Increment the last octet of the IP addresses to ensure unique subnets
  
  }
  

  // Second level dumbbell
  PointToPointDumbbellHelper secondLevelDumbbell(numberOfCentralNodes, p2p, 1, p2p, p2p);

  for (uint32_t i = 0; i < secondLevelDumbbell.LeftCount(); ++i)
  {
    stack.Install(secondLevelDumbbell.GetLeft(i)); 
  }
  for (uint32_t i = 0; i < secondLevelDumbbell.RightCount(); ++i)
  {
    stack.Install(secondLevelDumbbell.GetRight(i));
  }
  stack.Install(secondLevelDumbbell.GetLeft());
  stack.Install(secondLevelDumbbell.GetRight());

  secondLevelDumbbell.AssignIpv4Addresses(leftIp, rightIp, routerIp);
  
  leftIp.NewNetwork ();
  rightIp.NewNetwork ();
  routerIp.NewNetwork ();

  Ipv4StaticRoutingHelper staticRoutingHelper;

  for (uint32_t i = 0; i < numberOfCentralNodes; ++i) 
  {
    NetDeviceContainer devices = p2p.Install(secondLevelDumbbell.GetLeft(i), firstLevelDumbbells[i].GetRight(0));
    Ipv4InterfaceContainer internetIpIfaces = ipv4h.Assign(devices);

    ipv4h.NewNetwork();

    Ptr<Ipv4> ipv4_l = secondLevelDumbbell.GetLeft(i)->GetObject<Ipv4>();
    Ipv4Address addr_l = ipv4_l->GetAddress(1,0).GetLocal();

    Ptr<Ipv4> ipv4_r = firstLevelDumbbells[i].GetRight(0)->GetObject<Ipv4>();
    Ipv4Address addr_r = ipv4_r->GetAddress(1,0).GetLocal();

    /*

    Ptr<Ipv4StaticRouting> staticRouting_l = staticRoutingHelper.GetStaticRouting(secondLevelDumbbell.GetLeft(i)->GetObject<Ipv4>());
    Ptr<Ipv4StaticRouting> staticRouting_r = staticRoutingHelper.GetStaticRouting(firstLevelDumbbells[i].GetRight(0)->GetObject<Ipv4>());

    // Assuming you want to route to the networks associated with these IPs
    Ipv4Mask networkMask("255.255.255.0");
    staticRouting_l->PrintRoutingTable(Create<OutputStreamWrapper>(&std::cout));
    staticRouting_r->PrintRoutingTable(Create<OutputStreamWrapper>(&std::cout));

        // Add static routes as needed
    staticRouting_l->AddNetworkRouteTo(Ipv4Address(addr_r.Get() & networkMask.Get()), networkMask, devices.Get(0)->GetIfIndex());
    staticRouting_r->AddNetworkRouteTo(Ipv4Address(addr_l.Get() & networkMask.Get()), networkMask, devices.Get(1)->GetIfIndex());

    std::cout<<"Routed added: "<< devices.Get(0)->GetIfIndex() << " "<< devices.Get(1)->GetIfIndex() << std::endl;
    
    staticRouting_l->PrintRoutingTable(Create<OutputStreamWrapper>(&std::cout));
    staticRouting_r->PrintRoutingTable(Create<OutputStreamWrapper>(&std::cout));

    std::cout << "IP address of node is " << addr_l << " " << addr_r << std::endl;
    */
  }

  // Connect server node to the second level dumbbell right leaf
  NetDeviceContainer internetDevices = p2p.Install(serverNodes.Get(0), secondLevelDumbbell.GetRight(0));
  Ipv4InterfaceContainer internetIpIfaces = ipv4h.Assign(internetDevices);

  
  // Application installation
  uint16_t port = 8080;
  for (uint32_t i = 0 ; i < numberOfCentralNodes; i++) {
    for (uint32_t j = 0; j < numberOfLowerNodes; j++) {

        Ptr<Node> lowerNode = firstLevelDumbbells[i].GetLeft(j);
        Ptr<Node> upperNode = secondLevelDumbbell.GetLeft(i);//firstLevelDumbbells[i].GetRight(0);

        Ptr<Socket> recvSocket = Socket::CreateSocket (lowerNode,  TcpSocketFactory::GetTypeId ());
        Ptr<Socket> sendSocket = Socket::CreateSocket (lowerNode, TcpSocketFactory::GetTypeId ());
        Address sendAddress = InetSocketAddress (upperNode->GetObject<Ipv4>()->GetAddress(1,0).GetLocal(), port);

        NS_LOG_UNCOND("[First] Send address for Node "<<lowerNode->GetObject<Ipv4>()->GetAddress(1, 0).GetLocal()<<" is set as "<<InetSocketAddress::ConvertFrom(sendAddress).GetIpv4());

        Ptr<ReceiveSendApplication> app = CreateObject<ReceiveSendApplication> (ReceiveSendApplication::TRANSMIT_ONLY, 500, MilliSeconds(100), 0.1, 0.1);
        app->Setup (recvSocket, sendSocket, sendAddress, packetSize, interval);
        lowerNode->AddApplication (app);
        app->SetStartTime (Seconds (3.0 + i*0.1));
        app->SetStopTime (Seconds (4.0));
    }
  }

/*
  for (uint32_t i = 0; i < numberOfCentralNodes; i++) {
    Ptr<Node> lowerNode = firstLevelDumbbells[i].GetRight(0);//secondLevelDumbbell.GetLeft(i);
    Ptr<Node> upperNode = secondLevelDumbbell.GetLeft(i);

    Ptr<Socket> recvSocket = Socket::CreateSocket (lowerNode,  TcpSocketFactory::GetTypeId ());
    Ptr<Socket> sendSocket = Socket::CreateSocket (lowerNode, TcpSocketFactory::GetTypeId ());
    Address sendAddress = InetSocketAddress (upperNode->GetObject<Ipv4>()->GetAddress(1,0).GetLocal(), port);

    NS_LOG_UNCOND("[1.5] Send address for Node "<<lowerNode->GetObject<Ipv4>()->GetAddress(1, 0).GetLocal()<<" is set as "<<InetSocketAddress::ConvertFrom(sendAddress).GetIpv4());

    Ptr<ReceiveSendApplication> app = CreateObject<ReceiveSendApplication> (ReceiveSendApplication::RECEIVE_AND_TRANSMIT, 500,MilliSeconds(100), 0, 1);
    app->Setup (recvSocket, sendSocket, sendAddress, packetSize, interval);
    lowerNode->AddApplication (app);
    app->SetStartTime (Seconds (2.5 + i*0.1));
    app->SetStopTime (Seconds (4.0));
  }
*/

  for (uint32_t i = 0; i < numberOfCentralNodes; i++) {
    Ptr<Node> lowerNode = secondLevelDumbbell.GetLeft(i);
    Ptr<Node> upperNode = serverNodes.Get(0);

    Ptr<Socket> recvSocket = Socket::CreateSocket (lowerNode,  TcpSocketFactory::GetTypeId ());
    Ptr<Socket> sendSocket = Socket::CreateSocket (lowerNode, TcpSocketFactory::GetTypeId ());
    Address sendAddress = InetSocketAddress (upperNode->GetObject<Ipv4>()->GetAddress(1,0).GetLocal(), port);

    NS_LOG_UNCOND("[Second] Send address for Node "<<lowerNode->GetObject<Ipv4>()->GetAddress(1, 0).GetLocal()<<" is set as "<<InetSocketAddress::ConvertFrom(sendAddress).GetIpv4());

    Ptr<ReceiveSendApplication> app = CreateObject<ReceiveSendApplication> (ReceiveSendApplication::RECEIVE_AND_TRANSMIT, 500,MilliSeconds(100), 0.1, 0.1);
    app->Setup (recvSocket, sendSocket, sendAddress, packetSize, interval);
    lowerNode->AddApplication (app);
    app->SetStartTime (Seconds (2.0 + i*0.1));
    app->SetStopTime (Seconds (4.0));
  }
 
  Ptr<Node> lowerNode = serverNodes.Get(0); //secondLevelDumbbell.GetRight(0);

  Ptr<Socket> recvSocket = Socket::CreateSocket (lowerNode, TcpSocketFactory::GetTypeId ());
  Ptr<Socket> sendSocket = Socket::CreateSocket (lowerNode, TcpSocketFactory::GetTypeId ());

  NS_LOG_UNCOND("[Third] Send address for Node "<<lowerNode->GetObject<Ipv4>()->GetAddress(1, 0).GetLocal());

  Ptr<ReceiveSendApplication> app = CreateObject<ReceiveSendApplication> (ReceiveSendApplication::RECEIVE_ONLY, 500, MilliSeconds(0), 0.1, 0.1);
  app->Setup (recvSocket, sendSocket, InetSocketAddress (Ipv4Address::GetAny (), port), packetSize, interval); // Sending to an arbitrary address
  lowerNode->AddApplication (app);
  app->SetStartTime (Seconds (1.0));
  app->SetStopTime (Seconds (4.0));

  Ipv4GlobalRoutingHelper::PopulateRoutingTables();

  Ipv4GlobalRoutingHelper g;
  Ptr<OutputStreamWrapper> routingStream = Create<OutputStreamWrapper>("routes.txt", std::ios::out);
  g.PrintRoutingTableAllAt(Seconds(2), routingStream);
  
  std::string filename = "completed_objects.csv";
  Logger::Instance(filename).WriteCSVHeader();

  // Run the simulation
  Simulator::Run ();
  Simulator::Destroy ();

  return 0;
}