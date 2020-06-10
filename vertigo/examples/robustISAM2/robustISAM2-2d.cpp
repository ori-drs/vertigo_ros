/*
 * robustISAM2.cpp
 *
 *  Created on: Jul 13, 2012
 *      Author: niko
 */

#include "planarSLAM.h"
#include <gtsam/nonlinear/ISAM2.h>
#include <gtsam/nonlinear/NonlinearFactorGraph.h>
#include <gtsam/nonlinear/Marginals.h>

#include <iostream>
#include <fstream>

using namespace gtsam;

// additional boost features
#include "boost/program_options.hpp"
namespace po = boost::program_options;

#include "boost/foreach.hpp"
#define foreach BOOST_FOREACH

//#include <boost/progress.hpp>
#include <boost/lexical_cast.hpp>

#include "betweenFactorSwitchable.h"
#include "switchVariableLinear.h"
#include "switchVariableSigmoid.h"
#include "betweenFactorMaxMix.h"
#include "betweenFactorAdaptive.h"
#include "shapeParameter.h"
#include "priorFactorOutlierProcess.h"
#include "timer.h"
using namespace vertigo;


// ===================================================================
struct Pose {
  int id;
  double x,y,th;
};

struct Edge {
  int i, j;
  double x, y, th;
  bool switchable;
  bool maxMix;
  Matrix covariance;
  double weight;
};



void writeResults(Values &results, std::string outputFile)
{
  std::ofstream resultFile(outputFile.c_str());

  Values::ConstFiltered<Pose2> result_poses = results.filter<Pose2>();
  foreach (const Values::ConstFiltered<Pose2>::KeyValuePair& key_value, result_poses) {
    Pose2 p = key_value.value;
    resultFile << "VERTEX_SE2 " << p.x() << " " << p.y() << " " << p.theta() << endl;
  }

  {
    Values::ConstFiltered<SwitchVariableLinear> result_switches = results.filter<SwitchVariableLinear>();
    foreach (const Values::ConstFiltered<SwitchVariableLinear>::KeyValuePair& key_value, result_switches) {
      resultFile << "VERTEX_SWITCH " << key_value.value.value() << endl;
    }
  }

  {
    Values::ConstFiltered<SwitchVariableSigmoid> result_switches = results.filter<SwitchVariableSigmoid>();
    foreach (const Values::ConstFiltered<SwitchVariableSigmoid>::KeyValuePair& key_value, result_switches) {
      resultFile << "VERTEX_SWITCH " << key_value.value.value() << endl;
    }
  }
}

// ===================================================================
bool parseDataset(std::string inputFile, std::vector<Pose>&poses, std::vector<Edge>&edges,std::multimap<int, int> &poseToEdges)
{
   std::cout << endl << "Parsing dataset " << inputFile << " ... " << std::endl;

	 // open the file
   std::ifstream inFile(inputFile.c_str());
	 if (!inFile) {
     std::cerr << "Error opening dataset file " << inputFile << std::endl;
		 return false;
	 }

//   std::cout << "in file eof is: " << inFile.eof() << endl;
	 // go through the dataset file
	 while (!inFile.eof()) {
		 // depending on the type of vertex or edge, read the data
     std::string type;
		 inFile >> type;
//     std::cout << "in file is: " << type << endl;
		 if (type == "VERTEX_SE2") {
			 Pose p;
			 inFile  >> p.id >> p.x >> p.y >> p.th;
			 poses.push_back(p);
//       std::cout << "size of poses is: " << poses.size() << endl;
		 }

		 else if (type == "EDGE_SE2_SWITCHABLE" || type == "EDGE_SE2" || type == "EDGE_SE2_MAXMIX") {
			 double dummy;
			 Edge e;

			 // read pose IDs
			 inFile >> e.i >> e.j;

			 if (e.i>e.j) {
         std::swap(e.i,e.j);
			 }

			 // read the switch variable ID (only needed in g2o, we dont need it here in the gtsam world)
			 if (type == "EDGE_SE2_SWITCHABLE") inFile >> dummy;
			 if (type == "EDGE_SE2_MAXMIX") inFile >> e.weight;


			 // read odometry measurement
			 inFile >> e.x >> e.y >> e.th;

			 // read information matrix
			 double info[6];
			 inFile >> info[0] >> info[1] >> info[2] >> info[3] >> info[4] >> info[5];
       Eigen::Matrix3d informationMatrix; //= Matrix_(3,3, info[0], info[1], info[2], info[1], info[3], info[4], info[2], info[4], info[5]);
       informationMatrix(0,0) = info[0];
       informationMatrix(0,1) = info[1];
       informationMatrix(0,2) = info[2];
       informationMatrix(1,0) = info[1];
       informationMatrix(1,1) = info[3];
       informationMatrix(1,2) = info[4];
       informationMatrix(2,0) = info[2];
       informationMatrix(2,1) = info[4];
       informationMatrix(2,2) = info[5];

			 e.covariance = inverse(informationMatrix);

			 if (type == "EDGE_SE2_SWITCHABLE") {
			   e.switchable=true;
			   e.maxMix=false;
			 }
			 else if (type == "EDGE_SE2_MAXMIX") {
			   e.switchable=false;
			   e.maxMix=true;
			 }
			 else {
			   e.switchable=false;
			   e.maxMix=false;
			 }

			 edges.push_back(e);

			 int id = edges.size()-1;
       poseToEdges.insert(std::pair<int, int>(e.j, id));
		 }
		 // else just read the next item at next iteration until one of the if clauses is true
	 }

     return true;
}


// ===================================================================
int main(int argc, char *argv[])
{
    // === handle command line arguments ===
    std::string inputFile, outputFile;
    std::string method;
    int stop;
    double relinThresh;

    ISAM2Params isam2Params;


    po::options_description desc("Allowed options");
      desc.add_options()
      ("help,h", "Show this help message.")
      ("input,i", po::value<std::string>(&inputFile)->default_value("dataset.g2o"),"Load dataset from this file.")
      ("output,o", po::value<std::string>(&outputFile)->default_value("results.isam"),"Save results in this file.")
      ("stop", po::value<int>(&stop)->default_value(-1), "Stop after this many poses.")
      ("verbose,v", "verbose mode")
      ("linear", "Use the linear weight function.")
      ("sigmoid", "Use the sigmoid as the switch function Psi instead of a linear function.")
      ("adaptive", "Use the adaptive weight function.")
      ("dogleg", "Use Powell's dogleg method instead of Gauss-Newton.")
      ("qr", "Use QR factorization instead of Cholesky.")
      ("relinSkip", po::value<int>(&isam2Params.relinearizeSkip)->default_value(10), "Only relinearize any variables every relinearizeSkip calls to ISAM2::update (default: 10)")
      ("relinThresh", po::value<double>(&relinThresh)->default_value(0.1),"Only relinearize variables whose linear delta magnitude is greater than this threshold." )
      ("odoOnly", "Only use odometry edges, discard all other edges in the graph.")
    ;

    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv).options(desc).allow_unregistered().run(),vm);
    po::notify(vm);

    if (vm.count("help")) {
      cout << desc << "\n";
      return 1;
    }

    bool verbose;
    if (vm.count("verbose")) verbose=true;
    else verbose=false;

    bool useSigmoid;
    if (vm.count("sigmoid")) useSigmoid=true;
    else useSigmoid = false;

    bool useLinear;
    if (vm.count("linear")) useLinear=true;
    else useLinear = false;

    bool useAdaptive;
    if (vm.count("adaptive")) useAdaptive=true;
    else useAdaptive = false;

    // === read and parse input file ===
    std::vector<Pose> poses;
    std::vector<Edge> edges;
    std::multimap<int, int> poseToEdges;
    if (!parseDataset(inputFile, poses, edges, poseToEdges)) {
      std::cerr << "Error parsing the dataset."<< endl;
    	return 0;
    }
    cout << "Read " << poses.size() << " poses and " << edges.size() << " edges." << endl;


    // === set up iSAM2 ===
    //ISAM2Params isam2Params; already defined on top
    ISAM2DoglegParams doglegParams;
    ISAM2GaussNewtonParams gaussParams;

    if (vm.count("dogleg")) isam2Params.optimizationParams = doglegParams;
    else isam2Params.optimizationParams = gaussParams;

    if (vm.count("qr")) isam2Params.factorization = ISAM2Params::QR;

    if (vm.count("verbose")) {
    //  isam2Params.enableDetailedResults = true;
      isam2Params.evaluateNonlinearError = true;
    }
    isam2Params.relinearizeThreshold = relinThresh;


    ISAM2 isam2(isam2Params);


    // === go through the data and incrementally feed it into iSAM2 ===
    cout << "Running iSAM2 ..." << endl;

    // set up progress bar
    int progressLength = poses.size();
    if (stop>0 && stop<=poses.size()) progressLength=stop;
//    boost::progress_display show_progress(progressLength);

    int counter=0;
    int switchCounter=-1;

    planarSLAM::Values globalInitialEstimate;
    Timer timer;

    // iterate through the poses and incrementally build and solve the graph, using iSAM2
    foreach (Pose p, poses) {

      planarSLAM::Graph graph;
      planarSLAM::Values initialEstimate;


    	 // find all the edges that involve this pose

      timer.tic("findEdges");
     std::pair<std::multimap<int, int>::iterator, std::multimap<int, int>::iterator > ret = poseToEdges.equal_range(p.id);
      timer.toc("findEdges");

      for (std::multimap<int, int>::iterator it=ret.first; it!=ret.second; it++) {

    	  // look at the edge and see if it is switchable or not
    	  Edge e = edges[it->second];


    	  // see if this is an odometry edge, if yes, use it to initialize the current pose
    	  if (e.j==e.i+1) {
          timer.tic("initialize");
    	    Pose2 predecessorPose = isam2.calculateEstimate<Pose2>(planarSLAM::PoseKey(p.id-1));
    	    if (isnan(predecessorPose.x()) || isnan(predecessorPose.y()) || isnan(predecessorPose.theta())) {
    	      cout << "! Degenerated solution (NaN) detected. Solver failed." << endl;
    	      writeResults(globalInitialEstimate, outputFile);
            timer.print(cout);
    	      return 0;
    	    }
    	    initialEstimate.insertPose(p.id, predecessorPose * Pose2(e.x, e.y, e.th));
    	    globalInitialEstimate.insertPose(p.id,  predecessorPose * Pose2(e.x, e.y, e.th) );
          timer.toc("initialize");
    	  }

        timer.tic("addEdges");

    		if (!e.switchable && !e.maxMix) {
    		  if (!vm.count("odoOnly") || (e.j == e.i+1) ) {
    		    // this is easy, use the convenience functions of gtsam
    		    SharedNoiseModel odom_model = noiseModel::Gaussian::Covariance(e.covariance);
    		    graph.addOdometry(e.i, e.j, Pose2(e.x, e.y, e.th), odom_model);
    		  }
    		}
    		else if (e.switchable && !vm.count("odoOnly")) {
          if (useLinear) {
            // create new switch variable
            initialEstimate.insert(Symbol('s',++switchCounter),SwitchVariableLinear(1.0));

            // create switch prior factor
            SharedNoiseModel switchPriorModel = noiseModel::Diagonal::Sigmas(Vector1(1.0));
            boost::shared_ptr<PriorFactor<SwitchVariableLinear> > switchPriorFactor (new PriorFactor<SwitchVariableLinear> (Symbol('s',switchCounter), SwitchVariableLinear(1.0), switchPriorModel));
            graph.push_back(switchPriorFactor);


            // create switchable odometry factor
            SharedNoiseModel odom_model = noiseModel::Gaussian::Covariance(e.covariance);
            boost::shared_ptr<NonlinearFactor> switchableFactor(new BetweenFactorSwitchableLinear<Pose2>(planarSLAM::PoseKey(e.i), planarSLAM::PoseKey(e.j), Symbol('s', switchCounter), Pose2(e.x, e.y, e.th), odom_model));
            graph.push_back(switchableFactor);
    		  }
          else if (useSigmoid) {

    		    // create new switch variable
    		    initialEstimate.insert(Symbol('s',++switchCounter),SwitchVariableSigmoid(10.0));

    		    // create switch prior factor
            SharedNoiseModel switchPriorModel = noiseModel::Diagonal::Sigmas(Vector1(20.0));
    		    boost::shared_ptr<PriorFactor<SwitchVariableSigmoid> > switchPriorFactor (new PriorFactor<SwitchVariableSigmoid> (Symbol('s',switchCounter), SwitchVariableSigmoid(10.0), switchPriorModel));
    		    graph.push_back(switchPriorFactor);

    		    // create switchable odometry factor
    		    SharedNoiseModel odom_model = noiseModel::Gaussian::Covariance(e.covariance);
    		    boost::shared_ptr<NonlinearFactor> switchableFactor(new BetweenFactorSwitchableSigmoid<Pose2>(planarSLAM::PoseKey(e.i), planarSLAM::PoseKey(e.j), Symbol('s', switchCounter), Pose2(e.x, e.y, e.th), odom_model));
    		    graph.push_back(switchableFactor);
    		  }
          else if (useAdaptive) {
           switchCounter++;
            // create switch prior factor
            SharedNoiseModel adaptivePriorModel = noiseModel::Diagonal::Sigmas(Vector1(20.0));

            if (switchCounter == 0){
              initialEstimate.insert(planarSLAM::AlphaKey(), ShapeParameter(1));

            }


            // create switchable odometry factor
            SharedNoiseModel odom_model = noiseModel::Gaussian::Covariance(e.covariance);
            // robust error model
//            SharedNoiseModel odom_model_huber = noiseModel::Robust::Create(noiseModel::mEstimator::Huber::Create(1.345), odom_model);
            boost::shared_ptr<NonlinearFactor> adaptiveFactor(new BetweenFactorAdaptive<Pose2>(planarSLAM::PoseKey(e.i), planarSLAM::PoseKey(e.j), planarSLAM::AlphaKey(), Pose2(e.x, e.y, e.th), odom_model));
            graph.push_back(adaptiveFactor);

            if (switchCounter == 0)
            graph.add(PriorFactorOutlierProcess<Pose2>(planarSLAM::PoseKey(e.i), planarSLAM::PoseKey(e.j), planarSLAM::AlphaKey(), Pose2(e.x, e.y, e.th), odom_model, adaptivePriorModel));


          }
    		}
    		else if (e.maxMix && !vm.count("odoOnly")) {
    		  // create mixture odometry factor
    		  SharedNoiseModel odom_model = noiseModel::Gaussian::Covariance(e.covariance);
    		  SharedNoiseModel null_model = noiseModel::Gaussian::Covariance(e.covariance / e.weight);

    		  boost::shared_ptr<NonlinearFactor> maxMixFactor(new BetweenFactorMaxMix<Pose2>(planarSLAM::PoseKey(e.i), planarSLAM::PoseKey(e.j), Pose2(e.x, e.y, e.th), odom_model, null_model, e.weight));
    		  graph.push_back(maxMixFactor);
    		}
        timer.toc("addEdges");
    	}


    	if (p.id==0) {
    	  // add prior for first pose
        SharedDiagonal prior_model = noiseModel::Diagonal::Sigmas(Vector3( 0.01, 0.01, 0.01));
    		graph.addPrior(p.id, Pose2(p.x, p.y, p.th), prior_model);

    		// initial value for first pose
    		initialEstimate.insertPose(p.id, Pose2(p.x, p.y, p.th));
//        cout << "p is: " << p.x << " , " << p.y << " , " << p.th << endl;
    	}





    	if (verbose) {
        timer.tic("update");
    	  ISAM2Result result = isam2.update(graph, initialEstimate);
        timer.toc("update");
    	  cout << "cliques: " << result.cliques << "\terr_before: " << *(result.errorBefore) << "\terr_after: " << *(result.errorAfter) << "\trelinearized: " << result.variablesRelinearized << endl;
    	}
    	else  {
        timer.tic("update");
    	  isam2.update(graph, initialEstimate);
        cout << "counter: " << counter << endl;
        timer.toc("update");
    	}


//    	if (!verbose) ++show_progress;
    	if ( (counter++ >= stop) && (stop>0)) break;

    	if (false) {
        timer.tic("marginals");
    	  gtsam::Marginals marginals(isam2.getFactorsUnsafe(), isam2.getLinearizationPoint());
    	  gtsam::Matrix cov = marginals.marginalCovariance(gtsam::Symbol('x', p.id));
        timer.toc("marginals");
    	  cout << cov << endl << endl;
    	}

/*
    	if (counter % 100 == 0) {
    	  Values results = isam2.calculateEstimate();
    	  char fileName[1024];
    	  sprintf(fileName, "results/results_%05d.isam", counter);
    	  writeResults(results, string(fileName));
    	}
*/
    }


    // === do final estimation ===
    //cout << endl;
    Values results = isam2.calculateBestEstimate();

    results.print();


    timer.print(cout);

    // === write results ===

    cout << endl << endl << "Writing results." << endl;
    writeResults(results, outputFile);

}
