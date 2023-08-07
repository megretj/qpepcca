// inspired by sixth.cc in the tutorial section of ns3
// it doesn't work without the tutorial-app.h and tutorial-app.cc files in the same directory
// so put this file in the same folder as the tutorial files and then it is possible to run it.

#include "tutorial-app.h"

#include "ns3/applications-module.h"
#include "ns3/core-module.h"
#include "ns3/internet-module.h"
#include "ns3/network-module.h"
#include "ns3/point-to-point-module.h"

#include <fstream>
#include <iostream>
#include <cstdio>
#include <chrono>
#include <ctime>
#include <string>

using namespace ns3;

NS_LOG_COMPONENT_DEFINE("SixthScriptExample");

// ===========================================================================
//
//         node 0                 node 1
//   +----------------+    +----------------+
//   |    ns-3 TCP    |    |    ns-3 TCP    |
//   +----------------+    +----------------+
//   |    10.1.1.1    |    |    10.1.1.2    |
//   +----------------+    +----------------+
//   | point-to-point |    | point-to-point |
//   +----------------+    +----------------+
//           |                     |
//           +---------------------+
//                5 Mbps, 2 ms
//
//
// We want to look at changes in the ns-3 TCP congestion window.  We need
// to crank up a flow and hook the CongestionWindow attribute on the socket
// of the sender.  Normally one would use an on-off application to generate a
// flow, but this has a couple of problems.  First, the socket of the on-off
// application is not created until Application Start time, so we wouldn't be
// able to hook the socket (now) at configuration time.  Second, even if we
// could arrange a call after start time, the socket is not public so we
// couldn't get at it.
//
// So, we can cook up a simple version of the on-off application that does what
// we want.  On the plus side we don't need all of the complexity of the on-off
// application.  On the minus side, we don't have a helper, so we have to get
// a little more involved in the details, but this is trivial.
//
// So first, we create a socket and do the trace connect on it; then we pass
// this socket into the constructor of our simple application which we then
// install in the source node.
// ===========================================================================
//

/**
 * Congestion window change callback
 *
 * \param stream The output stream file.
 * \param oldCwnd Old congestion window.
 * \param newCwnd New congestion window.
 */
static void
CwndChange(Ptr<OutputStreamWrapper> stream, uint32_t oldCwnd, uint32_t newCwnd)
{
    // NS_LOG_UNCOND(Simulator::Now().GetSeconds() << "\t" << newCwnd);
    *stream->GetStream() << Simulator::Now().GetSeconds() << "\t" << newCwnd
                         << std::endl;
}
/**
 * Received packet callback
 * 
 * \param stream The output stream file.
 * \param p The received packet.
 * \param addr The address of the sender.
 */
static void
Rx(Ptr<OutputStreamWrapper> stream, Ptr<const Packet> p, const Address& addr)
{
    // NS_LOG_UNCOND(Simulator::Now().GetSeconds() << "\t" << p->GetSize());
    *stream->GetStream() << Simulator::Now().GetSeconds() << "\t" << p->GetSize()
                         << std::endl;
}

/**
 * Rx drop callback
 *
 * \param file The output PCAP file.
 * \param p The dropped packet.
 */
static void
RxDrop(Ptr<PcapFileWrapper> file, Ptr<const Packet> p)
{
    NS_LOG_UNCOND("RxDrop at " << Simulator::Now().GetSeconds());
    file->Write(Simulator::Now(), p);
}

static void
NextTxSequence(Ptr<OutputStreamWrapper> stream, SequenceNumber32 old [[maybe_unused]], SequenceNumber32 nextTx)
{
    NS_LOG_UNCOND("NextTxSequence at " << Simulator::Now().GetSeconds() << "\t" << nextTx);
}

static void
Measurement(std::string bandwidth = "1Mbps", std::string delay = "100ms", double error_p = 0.000001, std::string transport_prot = "TcpHybla", uint32_t iteration = 0, double beta= 0.7)
{
    std::string prefix_file_name = transport_prot+ "/" + transport_prot + std::to_string(error_p)+ "-" + delay+ "-" + bandwidth + "-" + std::to_string(iteration);
    
    // Config::SetDefault("ns3::TcpSocket::RcvBufSize", UintegerValue(1 << 22));
    // Config::SetDefault("ns3::TcpSocket::SndBufSize", UintegerValue(1 << 10));
    Config::SetDefault("ns3::TcpSocket::TcpNoDelay", BooleanValue(true));
    Config::SetDefault("ns3::TcpSocketBase::Sack", BooleanValue(false));
    Config::SetDefault("ns3::TcpHybla::RRTT", TimeValue(MilliSeconds(10)));
    Config::SetDefault("ns3::TcpCubic::Beta", DoubleValue(beta));
    Config::SetDefault("ns3::TcpCubic::FastConvergence", BooleanValue(false));
    Config::SetDefault("ns3::TcpCubic::HyStart", BooleanValue(false));
    Config::SetDefault("ns3::TcpL4Protocol::SocketType", StringValue("ns3::"+ transport_prot));
    Config::SetDefault("ns3::TcpSocket::InitialCwnd", UintegerValue(1));
    Config::SetDefault("ns3::TcpL4Protocol::RecoveryType",
                       TypeIdValue(TypeId::LookupByName("ns3::TcpClassicRecovery")));

    NodeContainer nodes;
    nodes.Create(2);

    PointToPointHelper pointToPoint;
    pointToPoint.SetDeviceAttribute("DataRate", StringValue(bandwidth));
    pointToPoint.SetChannelAttribute("Delay", StringValue(delay));

    NetDeviceContainer devices;
    devices = pointToPoint.Install(nodes);

    Ptr<RateErrorModel> em = CreateObject<RateErrorModel>();
    em->SetAttribute("ErrorRate", DoubleValue(error_p));
    devices.Get(1)->SetAttribute("ReceiveErrorModel", PointerValue(em));

    InternetStackHelper stack;
    stack.Install(nodes);

    Ipv4AddressHelper address;
    address.SetBase("10.1.1.0", "255.255.255.252");
    Ipv4InterfaceContainer interfaces = address.Assign(devices);

    uint16_t sinkPort = 8080;
    Address sinkAddress(InetSocketAddress(interfaces.GetAddress(1), sinkPort));
    PacketSinkHelper packetSinkHelper("ns3::TcpSocketFactory",
                                      InetSocketAddress(Ipv4Address::GetAny(), sinkPort));
    ApplicationContainer sinkApps = packetSinkHelper.Install(nodes.Get(1));
    sinkApps.Start(Seconds(0.));
    sinkApps.Stop(Seconds(100.));

    Ptr<Socket> ns3TcpSocket = Socket::CreateSocket(nodes.Get(0), TcpSocketFactory::GetTypeId());

    Ptr<TutorialApp> app = CreateObject<TutorialApp>();
    app->Setup(ns3TcpSocket, sinkAddress, 1040, 1000, DataRate("10Kbps")); // 1Gbps is the max rate from the sender
    nodes.Get(0)->AddApplication(app);
    app->SetStartTime(Seconds(1.));
    app->SetStopTime(Seconds(100.));

    AsciiTraceHelper asciiTraceHelper;
    Ptr<OutputStreamWrapper> stream = asciiTraceHelper.CreateFileStream("results/"+prefix_file_name+".cwnd");
    ns3TcpSocket->TraceConnectWithoutContext("CongestionWindow",
                                             MakeBoundCallback(&CwndChange, stream));
    ns3TcpSocket->TraceConnectWithoutContext("NextTxSequence", MakeBoundCallback(&NextTxSequence, stream));

    PcapHelper pcapHelper;
    Ptr<PcapFileWrapper> file =
        pcapHelper.CreateFile("sixth.pcap", std::ios::out, PcapHelper::DLT_PPP);
    devices.Get(1)->TraceConnectWithoutContext("PhyRxDrop", MakeBoundCallback(&RxDrop, file));

    Ptr<OutputStreamWrapper> stream2 = asciiTraceHelper.CreateFileStream("results/"+prefix_file_name+"-received.data");
    Config::ConnectWithoutContext("/NodeList/1/ApplicationList/*/$ns3::PacketSink/Rx", MakeBoundCallback(&Rx, stream2));

    Simulator::Stop(Seconds(100));
    Simulator::Run();
    Simulator::Destroy();
}

int
main(int argc, char* argv[])
{   
    CommandLine cmd(__FILE__);
    cmd.Parse(argc, argv);

    // Generate the logspace of error rates
    std::vector<double> errorRates = {0.000001};
    // double minError = std::pow(10, -5);
    // double maxError = std::pow(10, -0.01);
    // int numSteps = 30;
    // for (int i = 0; i < numSteps; ++i) {
    //     double err = minError * std::pow(maxError / minError, i / (numSteps - 1.0));
    //     errorRates.push_back(err);
    // }

    // Generate the linspace of delays = rtt/2 in ms
    // recall that the delays are doubled so the rtt is actually delay*2
    std::vector<double> delays;
    double minDelay = 5;
    double maxDelay = 300;
    int numSteps = 2;
    for (int i = 0; i < numSteps; ++i) {
        double delay = minDelay + (maxDelay - minDelay) * i / (numSteps - 1.0);
        delays.push_back(delay);
    }

    // Generate the linspace of bandwidths
    std::vector<std::string> bandwidths = {"1Mbps"};
    
    std::cout << "Starting simulation" << std::endl;
    auto start = std::chrono::system_clock::now();
    auto previous = start;
    // Iterate over the error rates and delays
    for (const auto& bandwidth :bandwidths){
        for (const auto& err : errorRates) {
            for (const auto& delay : delays) {
                for (int iter = 0; iter < 1; ++iter){
                // Measurement(bandwidth, std::to_string(delay) + "ms", err, "TcpHybla", iter);
                Measurement(bandwidth, std::to_string(delay) + "ms", err, "TcpCubic", iter);
                }
                std::cout << "Finished delay " << std::to_string(delay).c_str() << std::endl;
            }
            auto end = std::chrono::system_clock::now();
            std::chrono::duration<double> elapsed_seconds = end-start;
            std::chrono::duration<double> elapsed_since_previous = end-previous;
            std::cout << "Finished error rate " << std::to_string(err).c_str()
                    << " in " << elapsed_since_previous.count() << "s, "
                    << "elapsed time since begining: " << elapsed_seconds.count() << "s"
                    << std::endl;
            previous = end;
        }
        auto end = std::chrono::system_clock::now();
        std::chrono::duration<double> elapsed_seconds = end-start;
        std::cout  << "Finished bandwidth " << bandwidth.c_str()<< "elapsed time since begining: " << elapsed_seconds.count() << "s"
                    << std::endl;
    }

    return 0;
}