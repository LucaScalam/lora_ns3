/*
 * This script simulates a complex scenario with multiple gateways and end
 * devices. The metric of interest for this script is the throughput of the
 * network.
 */

#include "ns3/end-device-lora-phy.h"
#include "ns3/gateway-lora-phy.h"
#include "ns3/class-a-end-device-lorawan-mac.h"
#include "ns3/gateway-lorawan-mac.h"
#include "ns3/simulator.h"
#include "ns3/log.h"
#include "ns3/pointer.h"
#include "ns3/constant-position-mobility-model.h"
#include "ns3/lora-helper.h"
#include "ns3/node-container.h"
#include "ns3/mobility-helper.h"
#include "ns3/position-allocator.h"
#include "ns3/double.h"
#include "ns3/random-variable-stream.h"
#include "ns3/periodic-sender-helper.h"
#include "ns3/command-line.h"
#include "ns3/network-server-helper.h"
#include "ns3/correlated-shadowing-propagation-loss-model.h"
#include "ns3/building-penetration-loss.h"
#include "ns3/building-allocator.h"
#include "ns3/buildings-helper.h"
#include "ns3/forwarder-helper.h"
#include <algorithm>
#include <ctime>
#include <iostream>
#include <fstream>

#include "ns3/netanim-module.h"     //generate animation.xml file


using namespace ns3;
using namespace lorawan;

NS_LOG_COMPONENT_DEFINE ("ComplexLorawanNetworkExample");

// Network settings
int nDevices = 200;                     //Seteo #nodos_0
int N = 3;                              //Seteo #nodos_i a iterar
int nGateways = 1;
double radius = 2500;                   //Seteo radio_0
int R = 3;//4                              //Seteo #radio_j a iterar
double simulationTime = 1800;
int packSize = 50;                      //Seteo tamano de paquetes

// Channel model
bool realisticChannelModel = false;
// int appPeriodSeconds = 0;           //Aleatorio
int appPeriodSeconds = 600;

// Output control
bool print = true;


// CALCULO DE ToA
Time GetTimeOnAir (uint8_t sf, double bw, uint32_t pack_size, uint32_t nPreamble = 8)
{
  // Compute the symbol duration
  // Bandwidth is in Hz
  double tSym = pow (2, int(sf)) / (bw);

  // Compute the preamble duration
  double tPreamble = (double(nPreamble) + 4.25) * tSym;

  // Payload size
  //uint32_t pl = packet->GetSize ();      // Size in bytes
  //NS_LOG_DEBUG ("Packet of size " << pl << " bytes");

  // This step is needed since the formula deals with double values.
  // de = 1 when the low data rate optimization is enabled, 0 otherwise
  // h = 1 when header is implicit, 0 otherwise
  //double de = txParams.lowDataRateOptimizationEnabled ? 1 : 0;
  //double h = txParams.headerDisabled ? 1 : 0;
  //double crc = txParams.crcEnabled ? 1 : 0;

  // num and den refer to numerator and denominator of the time on air formula
  //double num = 8 * pl - 4 * txParams.sf + 28 + 16 * crc - 20 * h;
  double num = 8 * pack_size - 4 * sf + 28 + 16;
  //double den = 4 * (txParams.sf - 2 * de);
  double den = 4 * (sf);
  //double payloadSymbNb = 8 + std::max (std::ceil (num / den) *
  //                                     (txParams.codingRate + 4), double(0));
  double payloadSymbNb = 8 + std::max (std::ceil (num / den) * (1 + 4), double(0));

  // Time to transmit the payload
  double tPayload = payloadSymbNb * tSym;

  // Compute and return the total packet on-air time
  return Seconds (tPreamble + tPayload);
}



// MAIN
int
main (int argc, char *argv[])
{
  CommandLine cmd;
  cmd.AddValue ("nDevices", "Number of end devices to include in the simulation", nDevices);
  cmd.AddValue ("radius", "The radius of the area to simulate", radius);
  cmd.AddValue ("simulationTime", "The time for which to simulate", simulationTime);
  cmd.AddValue ("appPeriod",
                "The period in seconds to be used by periodically transmitting applications",
                appPeriodSeconds);
  cmd.AddValue ("print", "Whether or not to print various informations", print);
  cmd.Parse (argc, argv);

  // Set up logging
  // LogComponentEnable ("ComplexLorawanNetworkExample", LOG_LEVEL_ALL);
  // LogComponentEnable("LoraChannel", LOG_LEVEL_INFO);
  // LogComponentEnable("LoraPhy", LOG_LEVEL_ALL);
  // LogComponentEnable("EndDeviceLoraPhy", LOG_LEVEL_ALL);
  // LogComponentEnable("GatewayLoraPhy", LOG_LEVEL_ALL);
  // LogComponentEnable("LoraInterferenceHelper", LOG_LEVEL_ALL);
  // LogComponentEnable("LorawanMac", LOG_LEVEL_ALL);
  // LogComponentEnable("EndDeviceLorawanMac", LOG_LEVEL_ALL);
  // LogComponentEnable("ClassAEndDeviceLorawanMac", LOG_LEVEL_ALL);
  // LogComponentEnable("GatewayLorawanMac", LOG_LEVEL_ALL);
  // LogComponentEnable("LogicalLoraChannelHelper", LOG_LEVEL_ALL);
  // LogComponentEnable("LogicalLoraChannel", LOG_LEVEL_ALL);
  // LogComponentEnable("LoraHelper", LOG_LEVEL_ALL);
  // LogComponentEnable("LoraPhyHelper", LOG_LEVEL_ALL);
  // LogComponentEnable("LorawanMacHelper", LOG_LEVEL_ALL);
  // LogComponentEnable("PeriodicSenderHelper", LOG_LEVEL_ALL);
  // LogComponentEnable("PeriodicSender", LOG_LEVEL_ALL);
  // LogComponentEnable("LorawanMacHeader", LOG_LEVEL_ALL);
  // LogComponentEnable("LoraFrameHeader", LOG_LEVEL_ALL);
  // LogComponentEnable("NetworkScheduler", LOG_LEVEL_ALL);
  // LogComponentEnable("NetworkServer", LOG_LEVEL_ALL);
  // LogComponentEnable("NetworkStatus", LOG_LEVEL_ALL);
  // LogComponentEnable("NetworkController", LOG_LEVEL_ALL);

  



  /***********
   *  Setup  *
   ***********/

for (int i = 1; i < R+1; i++){                // itero sobre Radios
  for (int j = 1; j < N+1 ; j++){             // itero sobre #nodos
    //Abro archivo para guardar los resultados
    std::ofstream file;
    // file.open ("LoRa/Data/test_R"+std::to_string((int)(i*radius))+"_N"+std::to_string((int)(j*nDevices)) +".txt");    //Directorios Nahue
    file.open ("test_R"+std::to_string((int)(i*radius))+"_N"+std::to_string((int)(j*nDevices)) +".txt");       //Directorios Luca
    
    for (double k = 4; k < 7; k ++){        // itero sobre BW
      for (int l = 0; l<6; l++){            // itero sobre SF
        
  
  //Valor de BW en HZ
  int bw;
  switch((int)k)
  {
    case 4: 
      bw = 125000;
      break;
    case 5: 
      bw = 250000;
      break;
    case 6: 
      bw = 500000;
      break; 
  }

  //Forzamos el SF para obtener las métricas buscadas
  Config::SetDefault ("ns3::EndDeviceLorawanMac::DataRate", UintegerValue (l));        


  // Create the time value from the period
  Time appPeriod = Seconds (appPeriodSeconds);

  // Mobility
  MobilityHelper mobility;
  mobility.SetPositionAllocator ("ns3::UniformDiscPositionAllocator", "rho", DoubleValue (radius*i),
                                 "X", DoubleValue (0.0), "Y", DoubleValue (0.0));                     //SETEO RADIOS
  mobility.SetMobilityModel ("ns3::ConstantPositionMobilityModel");

  /************************
   *  Create the channel  *
   ************************/

  // Create the lora channel object
  Ptr<LogDistancePropagationLossModel> loss = CreateObject<LogDistancePropagationLossModel> ();
  loss->SetPathLossExponent (3.76);
  loss->SetReference (1, 7.7);

  if (realisticChannelModel)
    {
      // Create the correlated shadowing component
      Ptr<CorrelatedShadowingPropagationLossModel> shadowing =
          CreateObject<CorrelatedShadowingPropagationLossModel> ();

      // Aggregate shadowing to the logdistance loss
      loss->SetNext (shadowing);

      // Add the effect to the channel propagation loss
      Ptr<BuildingPenetrationLoss> buildingLoss = CreateObject<BuildingPenetrationLoss> ();

      shadowing->SetNext (buildingLoss);
    }

  Ptr<PropagationDelayModel> delay = CreateObject<ConstantSpeedPropagationDelayModel> ();

  Ptr<LoraChannel> channel = CreateObject<LoraChannel> (loss, delay);

  /************************
   *  Create the helpers  *
   ************************/

  // Create the LoraPhyHelper
  LoraPhyHelper phyHelper = LoraPhyHelper ();
  phyHelper.SetChannel (channel);

  // Create the LorawanMacHelper
  LorawanMacHelper macHelper = LorawanMacHelper ();
  // macHelper.Set (LorawanMacHelper::EU);                    //SETEO REGION
  //enum LorawanMacHelper::Regions {EU,SingleChannel,ALOHA};  //Regiones disponibles

  macHelper.SetRegion(LorawanMacHelper::EU);
  

  // Create the LoraHelper
  LoraHelper helper = LoraHelper ();
  helper.EnablePacketTracking (); // Output filename
  //helper.EnableSimulationTimePrinting ();

  //Create the NetworkServerHelper
  NetworkServerHelper nsHelper = NetworkServerHelper ();

  //Create the ForwarderHelper
  ForwarderHelper forHelper = ForwarderHelper ();

  /************************
   *  Create End Devices  *
   ************************/

  // Create a set of nodes
  NodeContainer endDevices;
  endDevices.Create (nDevices*j);                 //SETEO #NODOS

  // Assign a mobility model to each node
  mobility.Install (endDevices);

  // Make it so that nodes are at a certain height > 0
  for (NodeContainer::Iterator ji = endDevices.Begin (); ji != endDevices.End (); ++ji)
    {
      Ptr<MobilityModel> mobility = (*ji)->GetObject<MobilityModel> ();
      Vector position = mobility->GetPosition ();
      position.z = 1.2;
      mobility->SetPosition (position);
    }

  // Create the LoraNetDevices of the end devices
  uint8_t nwkId = 54;
  uint32_t nwkAddr = 1864;
  Ptr<LoraDeviceAddressGenerator> addrGen =
      CreateObject<LoraDeviceAddressGenerator> (nwkId, nwkAddr);

  // Create the LoraNetDevices of the end devices
  macHelper.SetAddressGenerator (addrGen);
  phyHelper.SetDeviceType (LoraPhyHelper::ED);
  macHelper.SetDeviceType (LorawanMacHelper::ED_A);
  helper.Install (phyHelper, macHelper, endDevices);


  // Now end devices are connected to the channel

  // Connect trace sources
  // for (NodeContainer::Iterator ji = endDevices.Begin (); ji != endDevices.End (); ++ji)
  //   {
  //     Ptr<Node> node = *ji;
  //     Ptr<LoraNetDevice> loraNetDevice = node->GetDevice (0)->GetObject<LoraNetDevice> ();
  //     Ptr<LoraPhy> phy = loraNetDevice->GetPhy ();
  //     // phy->SetChannel(channel);
  //   }

  /*********************
   *  Create Gateways  *
   *********************/

  // Create the gateway nodes (allocate them uniformely on the disc)
  NodeContainer gateways;
  gateways.Create (nGateways);

  Ptr<ListPositionAllocator> allocator = CreateObject<ListPositionAllocator> ();
  // Make it so that nodes are at a certain height > 0
  allocator->Add (Vector (0.0, 0.0, 15.0));
  mobility.SetPositionAllocator (allocator);
  mobility.Install (gateways);

  // Create a netdevice for each gateway
  phyHelper.SetDeviceType (LoraPhyHelper::GW);
  macHelper.SetDeviceType (LorawanMacHelper::GW);
  helper.Install (phyHelper, macHelper, gateways);

  /**********************
   *  Handle buildings  *
   **********************/

  double xLength = 130;
  double deltaX = 32;
  double yLength = 64;
  double deltaY = 17;
  int gridWidth = 2 * radius / (xLength + deltaX);
  int gridHeight = 2 * radius / (yLength + deltaY);
  if (realisticChannelModel == false)
    {
      gridWidth = 0;
      gridHeight = 0;
    }
  Ptr<GridBuildingAllocator> gridBuildingAllocator;
  gridBuildingAllocator = CreateObject<GridBuildingAllocator> ();
  gridBuildingAllocator->SetAttribute ("GridWidth", UintegerValue (gridWidth));
  gridBuildingAllocator->SetAttribute ("LengthX", DoubleValue (xLength));
  gridBuildingAllocator->SetAttribute ("LengthY", DoubleValue (yLength));
  gridBuildingAllocator->SetAttribute ("DeltaX", DoubleValue (deltaX));
  gridBuildingAllocator->SetAttribute ("DeltaY", DoubleValue (deltaY));
  gridBuildingAllocator->SetAttribute ("Height", DoubleValue (6));
  gridBuildingAllocator->SetBuildingAttribute ("NRoomsX", UintegerValue (2));
  gridBuildingAllocator->SetBuildingAttribute ("NRoomsY", UintegerValue (4));
  gridBuildingAllocator->SetBuildingAttribute ("NFloors", UintegerValue (2));
  gridBuildingAllocator->SetAttribute (
      "MinX", DoubleValue (-gridWidth * (xLength + deltaX) / 2 + deltaX / 2));
  gridBuildingAllocator->SetAttribute (
      "MinY", DoubleValue (-gridHeight * (yLength + deltaY) / 2 + deltaY / 2));
  BuildingContainer bContainer = gridBuildingAllocator->Create (gridWidth * gridHeight);

  BuildingsHelper::Install (endDevices);
  BuildingsHelper::Install (gateways);

  // Print the buildings
  if (print)
    {
      std::ofstream myfile;
      myfile.open ("buildings.txt");
      std::vector<Ptr<Building>>::const_iterator it;
      int ji = 1;
      for (it = bContainer.Begin (); it != bContainer.End (); ++it, ++ji)
        {
          Box boundaries = (*it)->GetBoundaries ();
          myfile << "set object " << ji << " rect from " << boundaries.xMin << "," << boundaries.yMin
                 << " to " << boundaries.xMax << "," << boundaries.yMax << std::endl;
        }
      myfile.close ();
    }


  /**********************************************
   *  Seteamos el BW mediante la mac de cada nodo  *
   **********************************************/
  
  std::vector<double> dataRate_V(7,bw);
  
  for (NodeContainer::Iterator ji = endDevices.Begin (); ji != endDevices.End (); ++ji)
    {
      Ptr<Node> node = *ji;
      Ptr<LoraNetDevice> loraNetDevice = node->GetDevice (0)->GetObject<LoraNetDevice> ();
      Ptr<LorawanMac> lmac = loraNetDevice->GetMac();
      lmac->SetBandwidthForDataRate (dataRate_V);
      // LorawanMacHelper::ApplyCommonEuConfigurations(lmac);
    }

  NS_LOG_DEBUG ("Completed configuration");


  /*********************************************
   *  Install applications on the end devices  *
   *********************************************/

  Time appStopTime = Seconds (simulationTime);
  PeriodicSenderHelper appHelper = PeriodicSenderHelper ();
  appHelper.SetPeriod (Seconds (appPeriodSeconds));
  appHelper.SetPacketSize (packSize);
  Ptr<RandomVariableStream> rv = CreateObjectWithAttributes<UniformRandomVariable> (
      "Min", DoubleValue (0), "Max", DoubleValue (10));
  ApplicationContainer appContainer = appHelper.Install (endDevices);

  appContainer.Start (Seconds (0));
  appContainer.Stop (appStopTime);

  /**************************
   *  Create Network Server  *
   ***************************/

  // Create the NS node
  NodeContainer networkServer;
  networkServer.Create (1);

  // Create a NS for the network
  nsHelper.SetEndDevices (endDevices);
  nsHelper.SetGateways (gateways);
  nsHelper.Install (networkServer);

  //Create a forwarder for each gateway
  forHelper.Install (gateways);

  ////////////////
  // Simulation //
  ////////////////
  Simulator::Stop (appStopTime + Hours (1));
  
  ///////////////////////////
  // NetAnim configuration //
  ///////////////////////////
 /*  AnimationInterface anim ("Lora.xml"); // Mandatory
  anim.EnablePacketMetadata (); */

  NS_LOG_INFO ("Running simulation...");
  Simulator::Run ();

  Simulator::Destroy ();

  ///////////////////////////
  // Print results to file //
  ///////////////////////////
  NS_LOG_INFO ("Computing performance metrics...");


  LoraPacketTracker &tracker = helper.GetPacketTracker ();


  ///////////////////////////////
  // Print results to terminal //
  ///////////////////////////////

  //Print
  std::cout <<"Radio: "+std::to_string((int)radius*i) 
            << "m, BW:" << (int)bw/1000 
            << "kHz, SF: "<< 12-l 
            << ", Nodes: "<< std::to_string(nDevices*j)
            << ", Enviados/Recibidos: " << tracker.CountMacPacketsGlobally (Seconds (0), appStopTime + Hours (1)) 
            << ", ToA: " << GetTimeOnAir(12-l,bw,packSize).GetSeconds()
            << "s" << std::endl;

  //Escribimos test_Rxxxx_Nxxx.txt
  file << std::to_string((int)radius*i) 
          << " " << bw 
          << " " << 12-l 
          << " " << std::to_string(nDevices*j)
          << " " << tracker.CountMacPacketsGlobally (Seconds (0), appStopTime + Hours (1)) 
          << " " << GetTimeOnAir(12-l,bw,packSize)
          << std::endl;
  
     }
    }
    file.close();
   }
  }
  
  return 0;
}