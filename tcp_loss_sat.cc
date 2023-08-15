#include "ns3/applications-module.h"
#include "ns3/core-module.h"
#include "ns3/enum.h"
#include "ns3/error-model.h"
#include "ns3/event-id.h"
#include "ns3/flow-monitor-helper.h"
#include "ns3/internet-module.h"
#include "ns3/ipv4-global-routing-helper.h"
#include "ns3/network-module.h"
#include "ns3/point-to-point-module.h"
#include "ns3/tcp-header.h"
#include "ns3/traffic-control-module.h"
#include "ns3/udp-header.h"

#include <iostream>
#include <cstdio>
#include <chrono>
#include <ctime>
#include <string>

using namespace ns3;

NS_LOG_COMPONENT_DEFINE("TcpComparison");

static std::map<uint32_t, bool> firstCwnd;                      //!< First congestion window.
static std::map<uint32_t, bool> firstSshThr;                    //!< First SlowStart threshold.
static std::map<uint32_t, Ptr<OutputStreamWrapper>> cWndStream; //!< Congstion window output stream.
static std::map<uint32_t, uint32_t> cWndValue;                      //!< congestion window value.
static std::map<uint32_t, uint32_t> ssThreshValue;                  //!< SlowStart threshold value.

/**
 * Congestion window tracer.
 *
 * \param context The context.
 * \param oldval Old value.
 * \param newval New value.
 */
static void
CwndTracer(std::string context, uint32_t oldval, uint32_t newval)
{
    uint32_t nodeId = 0;

    if (firstCwnd[nodeId])
    {
        *cWndStream[nodeId]->GetStream() << "0.0 " << oldval << std::endl;
        firstCwnd[nodeId] = false;
    }
    *cWndStream[nodeId]->GetStream() << Simulator::Now().GetSeconds() << " " << newval << std::endl;
    cWndValue[nodeId] = newval;
}

/**
 * Congestion window trace connection.
 *
 * \param cwnd_tr_file_name Congestion window trace file name.
 * \param nodeId Node ID.
 */
static void
TraceCwnd(std::string cwnd_tr_file_name, uint32_t nodeId)
{
    AsciiTraceHelper ascii;
    cWndStream[nodeId] = ascii.CreateFileStream(cwnd_tr_file_name);
    Config::Connect("/NodeList/" + std::to_string(nodeId) +
                        "/$ns3::TcpL4Protocol/SocketList/0/CongestionWindow",
                    MakeCallback(&CwndTracer));
}

static void
SetReceiveBuffer(uint32_t size, Ptr<Node> node)
{
    Config::Set("NodeList/0/$ns3::TcpSocket/RcvBufSize", UintegerValue(size));
    Config::Set("NodeList/0/$ns3::TcpSocket/SndBufSize", UintegerValue(size));
    NS_LOG_LOGIC("SetReceiveBuffer: " << size);
}

// static void
// Measurement(std::string bandwidth = "1Mbps", std::string delay = "100ms", double error_p = 0.01, std::string transport_prot = "TcpHybla", uint32_t iteration = 0, std::string access_bandwidth = "100Mbps", double beta= 0.7)
// {
int
main(int argc, char* argv[])
{   
    bool tracing = true;
    std::string prefix_file_name = "TcpVariantsSat";
    // std::string prefix_file_name = transport_prot+ "/" + transport_prot + std::to_string(error_p)+ "-" + delay+ "-" + bandwidth + "-" + std::to_string(iteration);
    uint64_t data_mbytes = 0; // 0 means unlimited sending capacity for BulkSendApplication
    uint32_t mtu_bytes = 400; // The MTU means the maximum transmission unit, which is the maximum size of a packet that can be transmitted over the network otherwise known as MSS (maximum segment size)
    double duration = 500.0;
    uint32_t run = 0;
    bool flow_monitor = false;
    bool sack = false;
    bool hystart = false; // Hybrid start
    std::string recovery = "ns3::TcpClassicRecovery";
    uint16_t port = 50000; // Port at which the receiver application (PacketSinkApplication on node 1) will listen for incoming packets
    std::string transport_prot = "TcpCubic";
    double error_p = 0.0;
    std::string bandwidth = "10Mbps";
    std::string delay = "0.100ms";
    uint64_t rtt0 = 25.0;
    uint16_t delAckCount = 1;
    bool pcap = false;
    uint32_t queuesize = 10;
    double beta = 0.7;
    uint32_t iteration = 0;
    bool debug = false;
    uint32_t sndBufSize= 1<<21;

    CommandLine cmd(__FILE__);
    cmd.AddValue("transport_prot",
                 "Transport protocol to use: TcpNewReno, TcpLinuxReno, "
                 "TcpHybla, TcpHighSpeed, TcpHtcp, TcpVegas, TcpScalable, TcpVeno, "
                 "TcpBic, TcpYeah, TcpIllinois, TcpWestwoodPlus, TcpLedbat, "
                 "TcpLp, TcpDctcp, TcpCubic, TcpBbr",
                 transport_prot);
    cmd.AddValue("error_p", "Packet error rate", error_p);
    cmd.AddValue("bandwidth", "Bottleneck bandwidth", bandwidth);
    cmd.AddValue("delay", "Bottleneck delay", delay);
    cmd.AddValue("tracing", "Flag to enable/disable tracing", tracing);
    cmd.AddValue("prefix_name", "Prefix of output trace file", prefix_file_name);
    cmd.AddValue("data", "Number of Megabytes of data to transmit", data_mbytes);
    cmd.AddValue("mtu", "Size of IP packets to send in bytes", mtu_bytes);
    cmd.AddValue("duration", "Time to allow flows to run in seconds", duration);
    cmd.AddValue("run", "Run index (for setting repeatable seeds)", run);
    cmd.AddValue("flow_monitor", "Enable flow monitor", flow_monitor);
    cmd.AddValue("pcap_tracing", "Enable or disable PCAP tracing", pcap);
    cmd.AddValue("sack", "Enable or disable SACK option", sack);
    cmd.AddValue("recovery", "Recovery algorithm type to use (e.g., ns3::TcpPrrRecovery", recovery);
    cmd.AddValue("rtt0", "Reference rtt for TcpHybla in ms", rtt0);
    cmd.AddValue("delAckCount", "Number of packets to wait before sending a TCP ACK", delAckCount);
    cmd.AddValue("queueSize", "Size of the queue in packets", queuesize);
    cmd.AddValue("debug","Enable debug output", debug);
    cmd.AddValue("SndBufSize", "Size of the send buffer in bytes", sndBufSize);
    cmd.Parse(argc, argv);
    transport_prot = std::string("ns3::") + transport_prot;

    //Setting the seed for nondeterministic results
    SeedManager::SetSeed(iteration+1);
    SeedManager::SetRun(run);

    // 2 MB of TCP buffer
    Config::SetDefault("ns3::TcpSocket::RcvBufSize", UintegerValue(1 << 21));
    Config::SetDefault("ns3::TcpSocket::SndBufSize", UintegerValue(1 << 21));
    Config::SetDefault("ns3::TcpSocket::TcpNoDelay", BooleanValue(true));
    Config::SetDefault("ns3::TcpSocketBase::Sack", BooleanValue(sack));
    Config::SetDefault("ns3::TcpHybla::RRTT", TimeValue(MilliSeconds(rtt0)));
    Config::SetDefault("ns3::TcpSocket::DelAckCount", UintegerValue(delAckCount));// Number of packets to wait before sending a TCP ACK

    Config::SetDefault("ns3::TcpCubic::Beta", DoubleValue(beta));
    Config::SetDefault("ns3::TcpCubic::FastConvergence", BooleanValue(false));
    Config::SetDefault("ns3::TcpCubic::HyStart", BooleanValue(hystart));

     Config::SetDefault("ns3::TcpL4Protocol::RecoveryType",
                       TypeIdValue(TypeId::LookupByName(recovery)));
    // Select TCP variant
    TypeId tcpTid;
    NS_ABORT_MSG_UNLESS(TypeId::LookupByNameFailSafe(transport_prot, &tcpTid),
                        "TypeId " << transport_prot << " not found");
    Config::SetDefault("ns3::TcpL4Protocol::SocketType",
                       TypeIdValue(TypeId::LookupByName(transport_prot)));

    NS_LOG_LOGIC("Selected TCP Variant");

    if(debug){
    // User may find it convenient to enable logging
        LogComponentEnable("TcpComparison", LOG_LEVEL_ALL);
        //LogComponentEnable("BulkSendApplication", LOG_LEVEL_INFO);
        //LogComponentEnable("Ipv4", LOG_LEVEL_ALL);
        //LogComponentEnable("Ipv4StaticRouting", LOG_LEVEL_ALL);
        //LogComponentEnable("PfifoFastQueueDisc", LOG_LEVEL_ALL);
        //LogComponentEnable("TcpVariantsComparison", LOG_LEVEL_ALL);
        //LogComponentEnable("TcpSocketBase", LOG_LEVEL_DEBUG);
        LogComponentEnable("TcpCubic", LOG_LEVEL_DEBUG);
        LogComponentEnable("TcpHybla", LOG_LEVEL_DEBUG);
        //LogComponentEnable("TcpCongestionOps", LOG_LEVEL_DEBUG);
    }

    // Calculate the ADU (application data unit) size
    Header* temp_header = new Ipv4Header();
    uint32_t ip_header = temp_header->GetSerializedSize();
    NS_LOG_LOGIC("IP Header size is: " << ip_header);
    delete temp_header;
    temp_header = new TcpHeader();
    uint32_t tcp_header = temp_header->GetSerializedSize();
    NS_LOG_LOGIC("TCP Header size is: " << tcp_header);
    delete temp_header;
    uint32_t tcp_adu_size = mtu_bytes - 20 - (ip_header + tcp_header);
    NS_LOG_LOGIC("TCP ADU size is: " << tcp_adu_size);

    // Set the simulation start and stop time
    double start_time = 0.1;
    double stop_time = start_time + duration;

    // Create gateways, sources, and sinks
    NodeContainer serverandclient;
    serverandclient.Create(2);

    // Configure the error model, here we use RateErrorModel with packet error rate
    Ptr<UniformRandomVariable> uv = CreateObject<UniformRandomVariable>();
    uv->SetStream(50);
    RateErrorModel error_model;
    error_model.SetRandomVariable(uv);
    error_model.SetUnit(RateErrorModel::ERROR_UNIT_PACKET);
    error_model.SetRate(error_p);

    //Install the internet stack on all the nodes
    InternetStackHelper stack;
    stack.InstallAll();

    uint32_t size = queuesize*mtu_bytes;
    // static_cast<uint32_t>((std::min(bottle_b).GetBitRate() / 8) *
    //                                       ((bottle_d) * 2).GetSeconds());

    Config::SetDefault("ns3::PfifoFastQueueDisc::MaxSize",
                       QueueSizeValue(QueueSize(QueueSizeUnit::PACKETS, size / mtu_bytes)));

    // Configure types of queue on the nodes
    TrafficControlHelper tchPfifo;
    tchPfifo.SetRootQueueDisc("ns3::PfifoFastQueueDisc");

    DataRate bottle_b(bandwidth);
    Time bottle_d(delay);

    // Configure the addresses of the nodes
    Ipv4AddressHelper address;
    address.SetBase("10.0.0.0", "255.255.255.0");

    // Configure a new channel between the source and the sink instead of the source-gateway and gateway sink
    PointToPointHelper SatelliteLink;
    SatelliteLink.SetDeviceAttribute("DataRate", StringValue(bandwidth));
    SatelliteLink.SetChannelAttribute("Delay", StringValue(delay));
    SatelliteLink.SetDeviceAttribute("ReceiveErrorModel", PointerValue(&error_model));

    NetDeviceContainer devices;
    devices = SatelliteLink.Install(serverandclient.Get(0), serverandclient.Get(1));
    tchPfifo.Install(devices);
    address.NewNetwork();
    Ipv4InterfaceContainer interfaces = address.Assign(devices);

    // Print the IP address of the sink and source nodes for debugging
    Ptr<Ipv4> ipv4 = serverandclient.Get(0)->GetObject<Ipv4>();
    Ipv4InterfaceAddress iaddr = ipv4->GetAddress(1, 0);
    Ipv4Address addri = iaddr.GetAddress();
    NS_LOG_INFO("The IP address of the source node is: " << addri << "or" << interfaces.GetAddress(0));
    ipv4 = serverandclient.Get(1)->GetObject<Ipv4>();
    iaddr = ipv4->GetAddress(1, 0);
    addri = iaddr.GetLocal();
    NS_LOG_INFO("The IP address of the sink node is: " << addri << " or " << interfaces.GetAddress(1));


    // Prepare Applications

    // First on the sender side, install a BulkSendApplication on node 0, which is the source
    AddressValue remoteAddress(InetSocketAddress(interfaces.GetAddress(1), port));
    Config::SetDefault("ns3::TcpSocket::SegmentSize", UintegerValue(tcp_adu_size));
    BulkSendHelper bulksender("ns3::TcpSocketFactory", Address());
    bulksender.SetAttribute("Remote", remoteAddress);
    bulksender.SetAttribute("SendSize", UintegerValue(tcp_adu_size));
    bulksender.SetAttribute("MaxBytes", UintegerValue(data_mbytes * 1000000));

    ApplicationContainer sourceApp = bulksender.Install(serverandclient.Get(0));
    sourceApp.Start(Seconds(0));
    sourceApp.Stop(Seconds(stop_time - 3));

    // Then on the receiver side, install a PacketSinkApplication on node 1, which is the sink
    Address sinkLocalAddress(InetSocketAddress(Ipv4Address::GetAny(), port));
    PacketSinkHelper sinkHelper("ns3::TcpSocketFactory", sinkLocalAddress);

    sinkHelper.SetAttribute("Protocol", TypeIdValue(TcpSocketFactory::GetTypeId()));
    ApplicationContainer sinkApp = sinkHelper.Install(serverandclient.Get(1));
    sinkApp.Start(Seconds(0));
    sinkApp.Stop(Seconds(stop_time));

    // Set up tracing if enabled
    if (tracing)
    {
        firstCwnd[0] = true;
        Simulator::Schedule(Seconds(0.00001),
                            &TraceCwnd,
                            "results/"+ prefix_file_name + "-cwnd.data",
                            0);
    }

    // Set up the buffer size of the receiver independently of the sender
    Simulator::Schedule(Seconds(0.00001),
                        &SetReceiveBuffer,
                        sndBufSize,
                        serverandclient.Get(1));

    // Flow monitor
    FlowMonitorHelper flowHelper;
    if (flow_monitor)
    {
        flowHelper.InstallAll();
    }
    // Now run the simulation
    Simulator::Stop(Seconds(stop_time));
    Simulator::Run();

    if (flow_monitor)
    {
        flowHelper.SerializeToXmlFile(prefix_file_name + ".flowmonitor", true, true);
    }

    Simulator::Destroy();
}

// int
// main(int argc, char* argv[])
// {   

//     // Generate the logspace of error rates
//     std::vector<double> errorRates = {0.000001, 0.01};
//     // double minError = std::pow(10, -5);
//     // double maxError = std::pow(10, -0.01);
//     // int numSteps = 30;
//     // for (int i = 0; i < numSteps; ++i) {
//     //     double err = minError * std::pow(maxError / minError, i / (numSteps - 1.0));
//     //     errorRates.push_back(err);
//     // }

//     // Generate the linspace of delays = rtt/2 in ms
//     // recall that the delays are doubled so the rtt is actually delay*2
//     std::vector<double> delays;
//     double minDelay = 5;
//     double maxDelay = 300;
//     int numSteps = 2;
//     for (int i = 0; i < numSteps; ++i) {
//         double delay = minDelay + (maxDelay - minDelay) * i / (numSteps - 1.0);
//         delays.push_back(delay);
//     }

//     // Generate the linspace of bandwidths
//     std::vector<std::string> bandwidths = {"1Mbps"};
    
//     std::cout << "Starting simulation" << std::endl;
//     auto start = std::chrono::system_clock::now();
//     auto previous = start;
//     // Iterate over the error rates and delays
//     for (const auto& bandwidth :bandwidths){
//         for (const auto& err : errorRates) {
//             for (const auto& delay : delays) {
//                 for (int iter = 0; iter < 1; ++iter){
//                 // Measurement(bandwidth, std::to_string(delay) + "ms", err, "TcpHybla", iter);
//                 Measurement(bandwidth, std::to_string(delay) + "ms", err, "TcpCubic",iter,"10Kbps");
//                 }
//                 std::cout << "Finished delay " << std::to_string(delay).c_str() << std::endl;
//             }
//             auto end = std::chrono::system_clock::now();
//             std::chrono::duration<double> elapsed_seconds = end-start;
//             std::chrono::duration<double> elapsed_since_previous = end-previous;
//             std::cout << "Finished error rate " << std::to_string(err).c_str()
//                     << " in " << elapsed_since_previous.count() << "s, "
//                     << "elapsed time since begining: " << elapsed_seconds.count() << "s"
//                     << std::endl;
//             previous = end;
//         }
//         auto end = std::chrono::system_clock::now();
//         std::chrono::duration<double> elapsed_seconds = end-start;
//         std::cout  << "Finished bandwidth " << bandwidth.c_str()<< "elapsed time since begining: " << elapsed_seconds.count() << "s"
//                     << std::endl;
//     }

//     return 0;
// }