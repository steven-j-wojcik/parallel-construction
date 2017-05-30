#include "graph.h"
#include <mpi.h>
#include <time.h>

#define master 0

int main(int argc, char * argv[]) 
{

	double trialsPP = 1000000;
	int myRank, numProcs; //trials per processor

	MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
    MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
    MPI_Status status;
    double totalTrials = trialsPP * (numProcs-1);
    //minus 1 because 1 processor does back tracking only. 
    //1 proc = serial, 2 procs = serial split between procs, 3+ procs = parallel



    srand(time(NULL)*myRank);	//needed for seed


    //srand(time(nullptr)*myRank);
    Graph G("input.txt");
    if (numProcs == 1 )
    {
    	double backstart = MPI_Wtime();
    	G.Backtrack(-1);
    	double backend = MPI_Wtime() - backstart;
		double constructStart = MPI_Wtime();
		G.Construction(trialsPP); 
		G.Calcval(trialsPP);
		double constructEnd = MPI_Wtime() - constructStart;
		G.setEstimatedRel();
		G.setExactMinCutPath();
		G.setEstimMinCutPath();
		G.Print_Num_PS_2(); //to be safe on pathsets, print both just incase something messes up in the 0.5 additiong thing.
		G.Print("dodec8np.txt");
		ofstream outFile;
   		outFile.open("dodec8np.txt", std::ios::app);
   		outFile << "=========================" << endl;
   		outFile << endl << "Number of Processors: " << numProcs << endl;
    	outFile << "Number of Trials per processor: " << trialsPP << endl;
    	outFile << "(proc 0 only does back tracking, doesnt count towards trials)" << endl;
    	outFile << "Total number of trials: " << trialsPP << endl;
    	outFile << "Backtracking time: " << backend << endl;
    	outFile << "Construction time: " << constructEnd << endl;
    	outFile << "Constructoin time is the time to aggragte all of the pathsets/totalTrials" << endl;
    	outFile.close();
    }
    else if (numProcs > 1)
    {
	    if (myRank == master)
	    {
	    	double backstart = MPI_Wtime();
			G.Backtrack(-1);
			double backend = MPI_Wtime() - backstart;
		    double * tempf = new double [G.getEdges()+1];
		    double constructStart = MPI_Wtime();
		    for(int i = 1; i < numProcs; i++)
		    {
		    	MPI_Recv(&tempf[0], G.getEdges(), MPI_DOUBLE, i, 0, MPI_COMM_WORLD, &status);
		    	for(int j = 0; j <  G.getEdges(); j++)
		    	{
		    		G.f[j] += tempf[j];
		    	}
		    }

		    G.Calcval(totalTrials);
		    G.setEstimatedRel();
			G.setExactMinCutPath();
			G.setEstimMinCutPath();
			double constructEnd = MPI_Wtime() - constructStart;
			G.Print_Num_PS_2(); //to be safe on pathsets, print both just incase something messes up in the 0.5 additiong thing.
		    G.Print("dodec4np.txt");
		    ofstream outFile;
		    outFile << "=========================" << endl;
    		outFile.open("dodec4np.txt", std::ios::app);
    		outFile << endl << "Number of Processors: " << numProcs << endl;
    		outFile << "Number of Trials per processor: " << trialsPP << endl;
    		outFile << "(proc 0 only does back tracking, doesnt count towards trials)" << endl;
    		outFile << "Total number of trials: " << totalTrials << endl;
    		outFile << "Backtracking time: " << backend << endl;
    		outFile << "Construction time: " << constructEnd << endl;
    		outFile << "Constructoin time is the time to aggragte all of the pathsets/totalTrials and get estim rel, mincut/path " << endl;
    		delete [] tempf;
    		outFile.close();
		}
		else
		{
	 		G.Construction(trialsPP);
			MPI_Send(&G.f[0], G.getEdges(), MPI_DOUBLE, master, 0, MPI_COMM_WORLD);
		}
	}
    MPI_Finalize();
	return 0;
} 