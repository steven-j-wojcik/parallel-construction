#include <iostream>
#include <fstream>
#include <iomanip>
#include <math.h>
#include <algorithm>
#include <limits>
using namespace std;

const int Max = 1500;

struct edge
{
	int v;
	int u;
	double prob;
};

class Graph
{
private:
	long double Rel; 			//overall reliability of graph (backtracking)
	long double EstimatedRel; 	//overall reliability of graph (from construction)
	int diameter;				//minimum path length between terminal nodes
	int nodes;					//number of nodes
	int edges;					//number of edges
	int edge_counter;			//temp varible used in backtracking
	edge * List; 				//set of edges
	int * Vector;				//temp varible
	int C[80][1500];
	int * Cut;					//unused
	int * Pathset;				//number of pathsets per subgraph size i accurate / backtracking
	double * f1;
	double * F;
	double * Cp;				//estimated pathsets per subgraph size i
	int * indices;				//temp varible used in construction
	int source;					//source node index
	int terminal;				//can be dyanamic array if you want more than 1 terminal.
	int counter;				//temp
	int counter2;				//temp
	double edge_rel;			//uniform edge reliability
	int minCut;
	int minPath;
	int estimatedMinCut;
	int estimatedMinPath;

public:
	double * f;
	~Graph()
	{
		delete [] List;
		delete [] Vector;
		delete [] Cut;
		delete [] Pathset;
		delete [] Cp;
		delete [] F;
		delete [] f;
		delete [] f1;
		delete [] indices;
	}
	Graph (Graph & other)
	{
		Rel = other.Rel;
		edge_rel = other.edge_rel;
		EstimatedRel = other.EstimatedRel;
		source = other.source;
		terminal = other.terminal;
		counter = 0;
		counter2 = 0;
		minCut = other.minCut;
		minPath = other.minPath;
		diameter = other.diameter;
		estimatedMinCut = other.estimatedMinCut;
		estimatedMinPath = other.estimatedMinPath;
		nodes = other.nodes;
		edges = other.edges;
		
		List=new edge[edges];
		Vector=new int[edges];
		Cut=new int[edges+1];
		Pathset=new int[edges+1];
		Cp = new double [edges+1];
		F = new  double [edges+1];
		f = new double [edges+1];
		f1 = new double [edges+1];
		indices= new int [edges+1];
		
		for (int j=0;j < edges; j++)
		{
			List[j].u = other.List[j].u;
			List[j].v = other.List[j].v;
			List[j].prob = other.List[j].prob;
			edge_rel = List[j].prob;
			Vector[j] = other.Vector[j];
			Cut[j] = other.Cut [j];
			Pathset[j] = other.Pathset[j];
			indices[j] = j;
			f[j]=other.f[j];
			f1[j]=other.f1[j];
			F[j]=other.F[j];
		}
		Cut[edges]=0;
		Pathset[edges]=0;
		indices[edges]=edges;
		f[edges]=0;
		f1[edges]=0;
		F[edges]=0;
	}
	Graph(int m, int D, double edgeRel)
	{
		if (D<=0)
		{
			nodes = 1;
			edges = 0;
			minCut = 0;
			minPath = 0;
			source = terminal = 0;
			return;
		}
		int bottomParrallel = 2*(D+1), topParrallel = D+1;
		Rel = 0.0;
		edge_rel = edgeRel;
		EstimatedRel = 0.0;
		source = 0;
		terminal = D;
		counter = 0;
		counter2 = 0;
		minCut = m;
		minPath = D;
		diameter = D;
		estimatedMinCut = 0;
		estimatedMinPath = 0;
		nodes = 3 * D + 4;
		edges = m * (3 * D + 9);;
		
		List=new edge[edges];
		Vector=new int[edges];
		Cut=new int[edges+1];
		Pathset=new int[edges+1];
		Cp = new double [edges+1];
		F = new  double [edges+1];
		f = new double [edges+1];
		f1 = new double [edges+1];
		indices= new int [edges+1];
		
		int edgeNum = 0;
		for (int i = 0; i < m; i++ )
		{
			//connect source to main body of gratph
			Create(edgeNum++, 0, 1, edgeRel);
			Create(edgeNum++, 0, topParrallel + 1, edgeRel);
			Create(edgeNum++, 0, bottomParrallel + 1, edgeRel);
			
			for(int j = 1; j <= D; j++)
			{
				Create(edgeNum++, j, j+1, edgeRel);
				Create(edgeNum++, j + topParrallel, j + topParrallel + 1, edgeRel);
				Create(edgeNum++, j + bottomParrallel, j + bottomParrallel + 1, edgeRel);
			}
			
			//vertical edges at front, if  D = 1, this connets the vertical edges of the terminal nodes
			Create(edgeNum++, 1, topParrallel + 1, edgeRel);
			Create(edgeNum++, 1, bottomParrallel + 1, edgeRel);
			
			//if the diameter isnt 1, add the vertical edges to the terminal node
			if (D != 1)
			{
				Create(edgeNum++, D, topParrallel + D, edgeRel);
				Create(edgeNum++, D, bottomParrallel + D, edgeRel);
			}
			//vertical edge at end
			Create(edgeNum++, D+1, topParrallel + D+1, edgeRel);
			Create(edgeNum++, D+1, bottomParrallel + D+1, edgeRel);
		}
		
		Cut[edges]=0;
		Pathset[edges]=0;
		indices[edges]=edges;
		f[edges]=0;
		f1[edges]=0;
		F[edges]=0;
	}
	Graph(int m, int Distance, int diam, double edgeRel)
	{
		if (Distance<=0)
		{
			nodes = 1;
			edges = 0;
			minCut = 0;
			minPath = 0;
			source = terminal = 0;
			return;
		}
		int bottomParrallel = 2*(Distance+1), topParrallel = Distance+1;
		Rel = 0.0;
		edge_rel = edgeRel;
		EstimatedRel = 0.0;
		source = 0;
		terminal = Distance;
		counter = 0;
		counter2 = 0;
		if (Distance == diam) minCut = m;
		else if (Distance > diam) minCut = 0;
		else if (Distance < diam) minCut = 3 * m;
		minPath = Distance;
		diameter = diam;
		estimatedMinCut = 0;
		estimatedMinPath = 0;
		nodes = 3 * Distance + 4;
		edges = m * (3 * Distance + 9);;
		
		List=new edge[edges];
		Vector=new int[edges];
		Cut=new int[edges+1];
		Pathset=new int[edges+1];
		Cp = new double [edges+1];
		F = new  double [edges+1];
		f = new double [edges+1];
		f1 = new double [edges+1];
		indices= new int [edges+1];
		
		int edgeNum = 0;
		for (int i = 0; i < m; i++ )
		{
			//connect source to main body of gratph
			Create(edgeNum++, 0, 1, edgeRel);
			Create(edgeNum++, 0, topParrallel + 1, edgeRel);
			Create(edgeNum++, 0, bottomParrallel + 1, edgeRel);
			
			for(int j = 1; j <= Distance; j++)
			{
				Create(edgeNum++, j, j+1, edgeRel);
				Create(edgeNum++, j + topParrallel, j + topParrallel + 1, edgeRel);
				Create(edgeNum++, j + bottomParrallel, j + bottomParrallel + 1, edgeRel);
			}
			
			//vertical edges at front, if  D = 1, this connets the vertical edges of the terminal nodes
			Create(edgeNum++, 1, topParrallel + 1, edgeRel);
			Create(edgeNum++, 1, bottomParrallel + 1, edgeRel);
			
			//if the diameter isnt 1, add the vertical edges to the terminal node
			if (Distance != 1)
			{
				Create(edgeNum++, Distance, topParrallel + Distance, edgeRel);
				Create(edgeNum++, Distance, bottomParrallel + Distance, edgeRel);
			}
			//vertical edge at end
			Create(edgeNum++, Distance + 1, topParrallel + Distance + 1, edgeRel);
			Create(edgeNum++, Distance + 1, bottomParrallel + Distance + 1, edgeRel);
		}
		
		Cut[edges]=0;
		Pathset[edges]=0;
		indices[edges]=edges;
		f[edges]=0;
		f1[edges]=0;
		F[edges]=0;
	}
	Graph (int n, int e,int d, int so, int te)
	{
		Rel = 0.0;
		EstimatedRel = 0.0;
		source = so;
		terminal = te;
		counter = 0;
		counter2 = 0;
		minCut = 0;
		minPath = 0;
		estimatedMinCut = 0;
		estimatedMinPath = 0;
		ifstream infile("input.txt");
		nodes = n;
		edges = e;
		diameter = d;
		List=new edge[edges];
		Vector=new int[edges];
		Cut=new int[edges+1];
		Pathset=new int[edges+1];
		Cp = new double [edges+1];
		F = new  double [edges+1];
		f = new double [edges+1];
		f1 = new double [edges+1];
		indices= new int [edges+1];
		
		for (int j=0;j < edges; j++)
		{
			infile >> List[j].u >> List[j].v;
			List[j].prob = 0.9;
			edge_rel = List[j].prob;
			Vector[j] = 1;
			Cut[j] = 0;
			Pathset[j] = 0;
			indices[j] = j;
			f[j]=0;
			f1[j]=0;
			F[j]=0;
		}
		infile.close();
		Cut[edges]=0;
		Pathset[edges]=0;
		indices[edges]=edges;
		f[edges]=0;
		f1[edges]=0;
		F[edges]=0;
	}
	
	Graph( char * fileName)
	{
		Rel = 0.0;
		EstimatedRel = 0.0;
		counter = 0;
		counter2 = 0;
		ifstream infile(fileName);
		infile >> diameter >> nodes >> edges;
		infile >> source >> terminal >> edge_rel;
		
		List=new edge[edges];
		Vector=new int[edges];
		Cut=new int[edges+1];
		Pathset=new int[edges+1];
		Cp = new double [edges+1];
		F = new  double [edges+1];
		f = new double [edges+1];
		f1 = new double [edges+1];
		indices= new int [edges+1];
		
		for (int j=0;j< edges; j++)
		{
			infile >> List[j].u >> List[j].v;
			List[j].prob = edge_rel;
			Vector[j] = 1;
			Cut[j] = 0;
			Pathset[j] = 0;
			indices[j] = j;
			f[j]=0.0;
			f1[j]=0.0;
			F[j]=0.0;
		}
		infile.close();
		Cut[edges]=0;
		Pathset[edges]=0;
		indices[edges]=edges;
		f[edges]=0.0;
		f1[edges]=0.0;
		F[edges]=0.0;
	}
	void Create();
	void Create(int edgeNum, int u, int v, double prob);
	void PrintCut (int j);
	void Update(int j, int Level);
	bool Diam(int d, int so, int te);
	void Backtrack(int edgenum);
	void Backtrack_mod(int edgenum, int Level);
	void PrintRel();
	void Print(string file);
	void Print_Num_PS (int s, int t);
	void Print_Num_PS ();
	void Print_Num_PS_2 ();
	void Print_Num_PS_2 (string file);
	void Print_Num_PS_3 ();
	bool Construction(double num_trials);
	long double Comb(int m, int l);
	bool Calcval (double num_trials);
	void SetST(int s1, int t1){source=s1; terminal= t1;Rel=0.0;}
	double getR(){return Rel;}
	int getEdges(){ return edges; }
	int getNodes(){ return nodes; }
	void setEstimatedRel();
	void setExactMinCutPath();
	void setEstimMinCutPath();
	int getEstimMinCut() { return estimatedMinCut; }
	int getEstimMinPath() { return estimatedMinPath; }
	int getExactMinCut() { return minCut; }
	int getExactMinPath() { return minPath; }
	void outputGraph()
	{
		cout << edges << ' ' << nodes << ' ' << diameter << ' ' << minCut << ' ' << minPath << endl;
		for(int i = 0; i < edges; i++)
		{
			cout<<List[i].u<<"->"<<List[i].v<<endl;
		}
	}
};

void Graph :: setEstimatedRel()
{
	for(int i = 1; i<=edges; i++)
	{
		EstimatedRel += Cp[i] * pow(edge_rel, i) * pow(1-edge_rel, edges-i);
	}
}

void Graph :: setExactMinCutPath()
{
	for( int i = edges ; i >= 1 ; i--)
	{
		if ( int(Comb(edges, i) - Pathset[i]) != 0)
		{
			minCut =  edges - i;
			break;
		}
	}
	for( int j = 1 ; j <= edges ; j++)
	{
		if (Pathset[j] != 0)
		{
			minPath = j;
			break;
		}
	}
}

void Graph :: setEstimMinCutPath()
{
	
	for( int i = edges ; i >= 1 ; i--)
	{
		if ( int(Comb(edges, i) - Cp[i]) != 0)
		{
			cout << i << ' ' << int(Comb(edges, i) - Cp[i]) <<endl;
			estimatedMinCut = edges - i;
			break;
		}
	}
	
	for( int j = 1 ; j <= edges ; j++)
	{
		if (Cp[j] != 0)
		{
			estimatedMinPath = j;
			break;
		}
	}
}

bool Graph:: Calcval (double num_trials)
{
	for (int i=1; i <= edges; i ++)
		f1[i] = double (f[i])/num_trials;
	int k = 1;
	double sumf = 0.0;
	for (int j = 1; j <= edges; j++)
	{
		sumf += f1[j];
		F[k] = sumf;
		Cp[k] = Comb(edges,k) * F[k];
		k++;
	}
	return (true);
}

long double Graph::Comb (int n, int j)
{
	double factorial = 1;
	if (j == 0)
		return (1);
	else
	{
		for (int i = n; i > (n-j) ; i--)
			factorial *=i ;
		for (int l = 1; l <= j; l ++)
			factorial /= l;
		return (factorial);
	}
}

void Graph::PrintCut (int j)
{
	cout << endl << "Cuts with " << j << "edges: " << Cut [j];
}

void Graph::Print_Num_PS (int s, int t)
{
	ofstream outfile;
	outfile.open("output.txt");		//rename with terminal, source, diameter, rel
	outfile << endl << "Ilya Graph";
	outfile << endl << "terminal vertices: " <<  s << ";" << t;
	outfile << endl << "diameter: " << diameter;
	outfile << endl << "Number of pathsets of given length ";
	outfile << endl;
	outfile << endl << "             0            1            2            3            4            5            6            7            8            9";
	outfile << endl << "-----------------------------------------------------------------------------------------------------------------------------------";
	for (int i=0; i <= (edges + 10); i = i + 10)
	{
		if (i != 0)
			outfile << endl << i/10 << "    ";
		else
			outfile << endl << "      ";// tenths
		for (int j=0; j < 10; j ++)
		{
			if ((i+j) <= edges)
				outfile << setw(8) << Pathset[i + j] << "     ";
			else
				break;
		}
	}
	outfile << endl;
	outfile << endl << "Number of pathsets calculated by contruction procedure";
	outfile << endl << "             0            1            2            3            4            5            6            7            8            9";
	outfile << endl << "-----------------------------------------------------------------------------------------------------------------------------------";
	for (int i=0; i <= (edges + 10); i = i + 10)
	{
		if (i != 0)
			outfile << endl << i/10 << "       ";
		else
			outfile << endl << "          ";// tenths
		for (int j=0; j < 10; j ++)
		{
			if ((i+j) <= edges)
				outfile << setw(8) << setprecision (2) << Cp[i + j]  << "     ";
			else
				break;
		}
	}
}

void Graph::Print_Num_PS ()
{
	ofstream outfile;
	outfile.open("output.txt");
	outfile << endl << "Ilya Graph";
	outfile << endl << "terminal vertices: " <<  source << ";" << terminal;
	outfile << endl << "diameter: " << diameter;
	outfile << endl << "Number of pathsets of given length ";
	outfile << endl;
	outfile << endl << "             0            1            2            3            4            5            6            7            8            9";
	outfile << endl << "-----------------------------------------------------------------------------------------------------------------------------------";
	for (int i=0; i <= (edges + 10); i = i + 10)
	{
		if (i != 0)
			outfile << endl << i/10 << "    ";
		else
			outfile << endl << "      ";// tenths
		for (int j=0; j < 10; j ++)
		{
			if ((i+j) <= edges)
				outfile << setw(8) << Pathset[i + j] << "     ";
			else
				break;
		}
	}
	outfile << endl;
	outfile << endl << "Number of pathsets calculated by contruction procedure";
	outfile << endl << "             0            1            2            3            4            5            6            7            8            9";
	outfile << endl << "-----------------------------------------------------------------------------------------------------------------------------------";
	for (int i=0; i <= (edges + 10); i = i + 10)
	{
		if (i != 0)
			outfile << endl << i/10 << "       ";
		else
			outfile << endl << "          ";// tenths
		for (int j=0; j < 10; j ++)
		{
			if ((i+j) <= edges)
				outfile << setw(8) << setprecision (2) << Cp[i + j]  << "     ";
			else
				break;
		}
	}
}

void Graph::Print_Num_PS_2 ()
{
	ofstream outfile;
	outfile.open("output.txt");
	outfile << endl << "Dodec";
	outfile << endl << "Edges: " << edges;
	outfile << endl << "Nodes: " << nodes;
	outfile << endl << "terminal vertices: " <<  source << ";" << terminal;
	outfile << endl << "diameter: " << diameter;
	outfile << endl << "uniform edge reliability: " << edge_rel;
	outfile << endl << "Number of pathsets of given length - Exact  Vs  Contruction";
	outfile << endl;
	outfile << endl << "Number of edges       # of Pathsets (Ex.)    # of Pathsets (Const)    % Diff.   ";
	outfile << endl << "---------------       -------------------    ---------------------    -------   ";
	
	for (int i=1; i <= edges; i++)
	{
		outfile << endl << setw (8) << i << "              " << setw(10) << Pathset[i] << "            " << setw(8) << fixed << setprecision(0) << (Cp[i]+0.5) << "            ";
		if (Pathset[i] > 0)
			outfile << setprecision(3) << setw(8) <<  float (int (Cp[i] + 0.5) - Pathset[i])/ Pathset[i] * 100.0;
		else
			outfile << setprecision(3) << setw(8) << 0.0;
	}
	outfile << endl << "-----------------------------------------------------------------" << endl;
	outfile << "\t\tmin path\tmin cut\treliability " << endl;
	outfile << "\t\t-------\t-------\t------------" << endl;
	outfile << "Exact:\t" << minPath << setw(11)<<' '<< minCut << setw(11) <<' '<< setprecision(20) << Rel << setprecision(1) <<endl;
	outfile << "Estim:\t" << estimatedMinPath << setw(11) <<' '<< estimatedMinCut << setw(11)<<' ' << setprecision(20) << EstimatedRel <<endl;
	outfile << "Precent Difference between exact and estimated reliability: " << setprecision(20) <<  ((Rel - EstimatedRel)/( (Rel + EstimatedRel)/2.0))*100.0 << " %" <<endl;
	outfile << "Precent Error between exact and estimated reliability: " << setprecision(20) << ((EstimatedRel - Rel)/Rel) * 100.0 << " %" <<endl;
	outfile.close();
}

void Graph::Print(string file)
{
	ofstream outfile;
	outfile.open(file.c_str());
	outfile << endl <<  file;
	outfile << endl << "Dodec";
	outfile << endl << "Edges: " << edges;
	outfile << endl << "Nodes: " << nodes;
	outfile << endl << "terminal vertices: " <<  source << ";" << terminal;
	outfile << endl << "diameter: " << diameter;
	outfile << endl << "uniform edge reliability: " << edge_rel;
	outfile << endl << "Number of pathsets of given length - Exact  Vs  Contruction";
	outfile << endl;
	outfile << endl << "Number of edges       # of Pathsets (Ex.)    # of Pathsets (Const)    % Diff.   ";
	outfile << endl << "---------------       -------------------    ---------------------    -------   ";
	
	for (int i=1; i <= edges; i++)
	{
		outfile << endl << setw (8) << i << "              " << setw(10) << fixed << setprecision(0) << Pathset[i]
		<< "            " << setw(8) << Cp[i] << "            ";
		if (Pathset[i] > 0)
			outfile << setprecision(3) << setw(8) <<  float (int (Cp[i] + 0.5) - Pathset[i])/ Pathset[i] * 100.0;
		else
			outfile << setprecision(3) << setw(8) << 0.0;
	}
	outfile << endl << "-----------------------------------------------------------------" << endl;
	outfile << "\t\tmin path\tmin cut\treliability " << endl;
	outfile << "\t\t-------\t-------\t------------" << endl;
	outfile << "Exact:\t" << minPath << setw(11)<<' '<< minCut << setw(11) <<' '<< setprecision(20) << Rel << setprecision(1) <<endl;
	outfile << "Estim:\t" << estimatedMinPath << setw(11) <<' '<< estimatedMinCut << setw(11)<<' ' << setprecision(20) << EstimatedRel <<endl;
	outfile << "Precent Difference between exact and estimated reliability: " << setprecision(20) <<  ((Rel - EstimatedRel)/( (Rel + EstimatedRel)/2.0))*100.0 << " %" <<endl;
	outfile << "Precent Error between exact and estimated reliability: " << setprecision(20) << ((EstimatedRel - Rel)/Rel) * 100.0 << " %" <<endl;;
	outfile.close();
}

void Graph::Print_Num_PS_2 (string file)
{
	ofstream outfile;
	outfile.open(file.c_str());
	outfile << endl <<  file;
	outfile << endl << "Edges: " << edges;
	outfile << endl << "Nodes: " << nodes;
	outfile << endl << "terminal vertices: " <<  source << ";" << terminal;
	outfile << endl << "diameter: " << diameter;
	outfile << endl << "uniform edge reliability: " << edge_rel;
	outfile << endl << "Number of pathsets of given length - Exact  Vs  Contruction";
	outfile << endl;
	outfile << endl << "Number of edges       # of Cutsets (Estim)    # of Pathsets (Estim)   ";
	outfile << endl << "---------------       -------------------     --------------------- ";
	
	for (int i=1; i <= edges; i++)
	{
		outfile << endl << setw (8) << i << "              " << setw(10) << fixed << setprecision(0) << Comb(edges, i) - Cp[i]
		<< "            " << setw(8) << Cp[i] << "            ";
	}
	outfile << endl << "-----------------------------------------------------------------" << endl;
	outfile << "\t\tmin path\tmin cut " << endl;
	outfile << "\t\t-------\t-------"<< endl;
	outfile << "Exact:\t" << minPath << "\t" << minCut << endl;
	outfile << "Estim:\t" << estimatedMinPath << "\t" << estimatedMinCut << endl;
	outfile.close();
}


void Graph::PrintRel()
{
	cout << setprecision(25);
	cout << endl << "The total Reliability with diameter < = " << diameter << " is " << Rel;
	cout << endl << "unique edge reliability: " << edge_rel;
}

bool Graph::Construction(double num_trials)
{
	int k=0;
	for (double numtrials=1; numtrials <= num_trials; numtrials ++)
	{
		for (int i=0; i < nodes; i++)
			for (int j=0 ; j < nodes; j++)
				if (i != j)
					C[i][j] = Max;
				else
					C[i][j]=0;
		random_shuffle(&indices[0], &indices[edges]);
		k=0;
		for (int l=0; l<edges; l++)
		{
			C[List[indices[l]].u][List[indices[l]].v]=1;
			C[List[indices[l]].v][List[indices[l]].u]=1; //  new edge up
			k++;
			if (this->Diam(diameter,source,terminal)) break; // add edges until the system is up
		}
		f[k]++;
	}
	return (true);
}

void Graph::Backtrack(int edgenum)
{
	long double Ri = 1.0;
	for (int i=0; i < nodes; i++)
	{
		for (int j=0 ; j < nodes; j++)
			if (i != j)
				C[i][j]=Max;
			else
				C[i][j]=0;
	}
	for (int i=0; i < edges; i++)
	{
		if (Vector[i])
		{
			C[List[i].u][List[i].v]=1;
			C[List[i].v][List[i].u]=1;
		}
	}
	edge_counter = 0;
	if (this->Diam(diameter,source,terminal))
	{
		counter ++; // new pathset
		for (int i=0; i < edges; i++)
			if (Vector[i])
			{
				Ri *= List[i].prob;
				edge_counter ++;
			}
			else
				Ri = Ri * (1.0 - List[i].prob);
		Rel += Ri; //add probability of occurance of operating state to the reliability
		Pathset[edge_counter]++;
		for (int j=edgenum+1; j < edges; j++)
		{
			Vector[j]=0;
			this->Backtrack(j);
			Vector[j]=1;
		}
	}
	else
		counter2 ++; // new cutset
}

//floyd warshall to check if the terminal is within diamter, D
bool Graph::Diam(int d, int so, int te)
{
	for (int k = 0 ; k < nodes; k++)
		for (int i = 0; i < nodes; i ++)
			for (int j = 0 ; j < nodes; j++)
				if((C[i][k] + C[k][j]) < C[i][j])
					C[i][j] = C[i][k] + C[k][j];
	if (C[so][te] > d)
		return (false);
	return (true);
}

void Graph::Create()
{
	for (int j=0; j < edges; j++)
	{
		cout << "Enter endpoint u: ";
		cin >> List[j].u;
		cout << endl;
		cout << "Enter endpoint v: ";
		cin >> List[j].v ;
		cout << endl;
	}
	cout << endl;
}

void Graph::Create(int edgeNum, int u, int v, double prob)
{
	Vector[edgeNum] = 1;
	Cut[edgeNum] = 0;
	Pathset[edgeNum] = 0;
	indices[edgeNum] = edgeNum;
	f[edgeNum]=0;
	f1[edgeNum]=0;
	F[edgeNum]=0;
	List[edgeNum].prob = prob;
	List[edgeNum].u = u;
	List[edgeNum].v = v;
}
