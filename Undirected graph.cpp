#include<iostream>
#include<vector>
#include<queue>
#include<new>
#include<algorithm>
#include<ctime>
#include<cstdlib>


using namespace std;


typedef long long LL;

typedef vector <LL> VLL;

typedef pair <LL,LL> PLL;

typedef vector <PLL> VPLL;

typedef priority_queue < PLL > PQPLL;


#define REP(i,a,b) for(int i=a; i<=b; ++i)


//const LL INF = (LL)1000000099 * 1000000009;

//const int UNDEF = -1;





class Vertex{
	private:
		int vertexID;
		
	public:
		Vertex(int vertexID = 0){
			this -> vertexID = vertexID;
		}
		
		~Vertex() {}
		
		Vertex(const Vertex &other){
			this -> vertexID = other.vertexID;
		}
		
		int getVertexID(){
			return vertexID;
		}
		
		void setVertexID(int vertexID){
			this -> vertexID = vertexID;
		}
		
		static Vertex createVertex(){
			int vertexID;
			cin >> vertexID;
			
			Vertex vertex(vertexID);
			return vertex;
		}
		
		static Vertex createRandomVertex(int maxVertexID){
			
			int vertexID = 1 + rand() % maxVertexID;
			
			Vertex vertex(vertexID);
			return vertex;
		}
};


class DirectedEdge{
	private:
		int directedEdgeID;
		Vertex startVertex;
		Vertex endVertex;
		int weight;
		
	public:
		DirectedEdge(){}
		
		
		DirectedEdge(int directedEdgeID, Vertex startVertex, Vertex endVertex, int weight = 1){
			this -> directedEdgeID = directedEdgeID;
			this -> startVertex = startVertex;
			this -> endVertex = endVertex;
			this -> weight = weight;
		}
		
		
		~DirectedEdge() {}
		
		
		DirectedEdge(const DirectedEdge &other){
			this -> directedEdgeID = other.directedEdgeID ;
			this -> startVertex = other.startVertex;
			this -> endVertex = other.endVertex;
			this -> weight = other.weight;
		}
		
		
		int getDirectedEdgeID(){
			return directedEdgeID;
		}
		
		
		void setDirectedEdgeID(int directedEdgeID){
			this -> directedEdgeID = directedEdgeID;
		}
		
		
		Vertex getStartVertex(){
			return startVertex;
		}
		
		
		void setStartVertex(Vertex startVertex){
			this -> startVertex = startVertex;
		}
		
		
		Vertex getEndVertex(){
			return endVertex;
		}
		
		
		void setEndVertex(Vertex endVertex){
			this -> endVertex = endVertex;
		}
		
		
		int getWeight(){
			return weight;
		}
		
		
		void setWeight(int weight){
			this -> weight = weight;
		}
		
		
		static DirectedEdge createWeightedDirectedEdge(int directedEdgeID){
			Vertex startVertex = Vertex::createVertex();
			Vertex endVertex = Vertex::createVertex();
			
			int weight;
			cin >> weight;
			
			DirectedEdge edge(directedEdgeID, startVertex, endVertex, weight);
			
			return edge;
		}
		
		
		static DirectedEdge createUnweightedDirectedEdge(int directedEdgeID){
			Vertex startVertex = Vertex::createVertex();
			Vertex endVertex = Vertex::createVertex();
			
			DirectedEdge edge(directedEdgeID, startVertex, endVertex);
			
			return edge;
		}
		
		
		void output(){
			int directedEdgeID = this -> directedEdgeID;
			int startVertex = (*this).getStartVertex().getVertexID();
			int endVertex = (*this).getEndVertex().getVertexID();
			int weight = this -> weight;
			
			cout << /*endl << directedEdgeID << " " <<*/ startVertex << " " << endVertex << " " << weight << endl;
		}
};


class UndirectedEdge{
	private:
		int undirectedEdgeID;
		DirectedEdge directedEdge1;
		DirectedEdge directedEdge2;
		
	public:
		UndirectedEdge() {}
		
		
		UndirectedEdge(int undirectedEdgeID, DirectedEdge directedEdge1, DirectedEdge directedEdge2){
			this -> undirectedEdgeID = undirectedEdgeID;
			this -> directedEdge1 = directedEdge1;
			this -> directedEdge2 = directedEdge2;
		}
		
		
		~UndirectedEdge() {}
		
		
		UndirectedEdge(const UndirectedEdge &other){
			this -> undirectedEdgeID = other.undirectedEdgeID;
			this -> directedEdge1 = other.directedEdge1;
			this -> directedEdge2 = other.directedEdge2;
		}
		
		
		int getUndirectedEdgeID(){
			return undirectedEdgeID;
		}
		
		
		void setUndirectedEdgeID(int undirectedEdgeID){
			this -> undirectedEdgeID = undirectedEdgeID;
		}
		
		
		DirectedEdge getDirectedEdge1(){
			return directedEdge1;
		}
		
		
		void setDirectedEdge1(DirectedEdge directedEdge1){
			this -> directedEdge1 = directedEdge1;
		}
		
		
		DirectedEdge getDirectedEdge2(){
			return directedEdge2;
		}
		
		
		void setDirectedEdge2(DirectedEdge directedEdge2){
			this -> directedEdge2 = directedEdge2;
		}
		
		
		int getWeight(){
			return directedEdge1.getWeight();
		}
		
		
		void setWeight(int weight){
			directedEdge1.setWeight(weight);
			directedEdge2.setWeight(weight);
		}
		
		
		static UndirectedEdge createWeightedUndirectedEdge(int undirectedEdgeID){
			Vertex startVertex = Vertex::createVertex();
			Vertex endVertex = Vertex::createVertex();
			
			int weight;
			cin >> weight;
			
			DirectedEdge directedEdge1( undirectedEdgeID,  startVertex,    endVertex,  weight );
			DirectedEdge directedEdge2( undirectedEdgeID,    endVertex,  startVertex,  weight );
			
			UndirectedEdge edge(undirectedEdgeID, directedEdge1, directedEdge2);
			
			return edge;
		}
		
		
		static UndirectedEdge createRandomWeightedUndirectedEdge(int undirectedEdgeID, int maxVertexID, int maxWeight){
			
			Vertex startVertex = Vertex::createRandomVertex(maxVertexID);
			Vertex endVertex = Vertex::createRandomVertex(maxVertexID);
			
			int weight = rand() % maxWeight;
			
			DirectedEdge directedEdge1( undirectedEdgeID,  startVertex,    endVertex,  weight );
			DirectedEdge directedEdge2( undirectedEdgeID,    endVertex,  startVertex,  weight );
			
			UndirectedEdge edge(undirectedEdgeID, directedEdge1, directedEdge2);
			
			return edge;
		}
		
		
		static UndirectedEdge createUnweightedUndirectedEdge(int undirectedEdgeID){
			Vertex startVertex = Vertex::createVertex();
			Vertex endVertex = Vertex::createVertex();
			
			DirectedEdge directedEdge1( 2*undirectedEdgeID - 1,  startVertex,    endVertex );
			DirectedEdge directedEdge2(     2*undirectedEdgeID,    endVertex,  startVertex );
			
			UndirectedEdge edge(undirectedEdgeID, directedEdge1, directedEdge2);
			
			return edge;
		}
		
		
		void output(){
			//int undirectedEdgeID = this -> undirectedEdgeID;
			DirectedEdge directedEdge1 =  this -> directedEdge1;
			//DirectedEdge directedEdge2 =  this -> directedEdge2;
			
			//cout << undirectedEdgeID << " " << startVertex << " " << endVertex << " " << weight << endl;
			directedEdge1.output();
			//directedEdge2.output();
			//cout << endl;
		}
};


class FindUnion{
	private:
		int N;
		int *repr;
		int *howMany;
		
	public:
		FindUnion(int N){
			this -> N = N;
			this -> repr = new int [N];
			this -> howMany = new int [N];
			
			REP(i,0,N-1){
				repr[i] = i;
				howMany[i] = 1;
			}
		}
		
		~FindUnion(){
			delete [] repr;
			delete [] howMany;
		}
		
		int Find(int x){
			int res = x;
			
			while( repr[res] != res ){
				res = repr[res];
			}
			
			repr[x] = res;
			
			return res;
		}
		
		void Union(int x, int y){
			int fx = Find(x), fy = Find(y);
			
			if( howMany[fx] > howMany[fy] ){
				repr[fy] = fx;
				howMany[fx] += howMany[fy];
			}
			else{
				repr[fx] = fy;
				howMany[fy] += howMany[fx];
			}
		}
		
};


class WeightedUndirectedGraph
{
	private:
		int numberOfVertices;
		int numberOfEdges;
		vector <UndirectedEdge> Edges;
		vector <DirectedEdge> *Graph;
		
		
		void operator += (DirectedEdge e){
			Graph[ e.getStartVertex().getVertexID() - 1 ].push_back( e );
		}
		
		
		void operator += (UndirectedEdge e){
			Edges.push_back(e);
			(*this) += e.getDirectedEdge1();
			(*this) += e.getDirectedEdge2();
		}
		
		
		void DFS_Depth(vector <int> &depth, int u, int parent){
			int numberOfNeighbours = Graph[u-1].size(), v;
			
			REP(i,0,numberOfNeighbours-1){
				v = Graph[u-1][i].getEndVertex().getVertexID();
				
				if( v != parent ){
					DFS_Depth( depth, v, u );
					
					depth[u-1] = max( depth[u-1], 1 + depth[v-1] );
				}
			}	
		}
	
		
		void DFS_SubtreeSize(vector <int> &subtreeSize, int u, int parent){
			int numberOfNeighbours = Graph[u-1].size(), v;
			
			REP(i,0,numberOfNeighbours-1){
				v = Graph[u-1][i].getEndVertex().getVertexID();
				
				if( v != parent ){
					DFS_SubtreeSize( subtreeSize, v, u );
					
					subtreeSize[u-1] += subtreeSize[v-1];	
				}
			}
		}
		
	
		void BFS_Coherentness(vector <bool> &visited, int startVertex){
			queue <int> Q;
				
			int u, numberOfNeighbours, v;
			
			Q.push( startVertex );
			while( !Q.empty() ){
				
				u = Q.front();    --u;
				Q.pop();
				
				visited[u] = true;
				
				numberOfNeighbours = Graph[u].size();
					
				REP(i,0,numberOfNeighbours-1){
					v = Graph[u][i].getEndVertex().getVertexID();    --v;
					
					if( !visited[v] ){
						Q.push(v+1);
					}
				}
			}
		}
		
		/*static LL min(LL a, LL b){
			return (a<b) ? a : b ;
		}*/
	
	
	public:
		
		WeightedUndirectedGraph(){
			numberOfVertices = 0;
			numberOfEdges = 0;
			Graph = NULL;
		}
		
		
		~WeightedUndirectedGraph(){
			Edges.clear();
			
			REP(u,0,numberOfVertices-1){
				Graph[u].clear();
			}
			delete [] Graph;
		}
		
		
		int getNumberOfVertices(){
			return numberOfVertices;
		}
		
		
		int getNumberOfEdges(){
			return numberOfEdges;
		}
		
		
		void createWeightedUndirectedGraph(){
			
			int numberOfVertices;
			cin >> numberOfVertices;
			
			int numberOfEdges;
			cin >> numberOfEdges;
			
			this -> numberOfVertices = numberOfVertices;
			this -> numberOfEdges = numberOfEdges;
			
			//cout << this -> numberOfVertices << " " << this -> numberOfEdges << endl;
			
			Graph = new vector <DirectedEdge> [this -> numberOfVertices];
			
			REP(edgeID,1,this -> numberOfEdges){
				UndirectedEdge e = UndirectedEdge::createWeightedUndirectedEdge(edgeID);
				//cout << e.getStartVertex().getVertexID() << " " << e.getEndVertex().getVertexID() << " " << e.getWeight() << endl;
				*this += e;
			}
		}
		
		
		void createRandomWeightedUndirectedGraph(int maxNumberOfVertices, int maxNumberOfEdges, int maxWeight){
			
			//map < pair<int,int>, int > myMap;
			 			
			this -> numberOfVertices = 1 + rand() % maxNumberOfVertices;
			this -> numberOfEdges = rand() % maxNumberOfEdges;
			
			//cout << this -> numberOfVertices << " " << this -> numberOfEdges << endl;
			
			Graph = new vector <DirectedEdge> [this -> numberOfVertices];
			
			//int u, v;
			REP(edgeID,1,this -> numberOfEdges){
				UndirectedEdge e = UndirectedEdge::createRandomWeightedUndirectedEdge(edgeID, numberOfVertices, maxWeight);
				//cout << e.getStartVertex().getVertexID() << " " << e.getEndVertex().getVertexID() << " " << e.getWeight() << endl;
				
				/*u = e.getDirectedEdge1().getStartVertex().getVertexID();
				v = e.getDirectedEdge1().getEndVertex().getVertexID();
				
				if( myMap[ make_pair(u,v) ] == 0  &&  myMap[ make_pair(v,u) ] == 0 ){
					
					myMap[ make_pair(u,v) ] = myMap[ make_pair(v,u) ] = 1;
				*/	
					*this += e;	
				//}
			}
		}
		
		
		bool isCoherent(){
			
			vector <bool> visited(numberOfVertices);
			REP(u,0,numberOfVertices-1){
				visited[u] = false;
			}
			
			int startVertex = 1;
			BFS_Coherentness(visited, startVertex);
			
			REP(u,0,numberOfVertices-1){
				if( !visited[u] ){
					return false;
				}
			}
			
			return true;
		}
		
		
		vector <int> computeDepths(int root){
			
			int UNDEF = -1;
			
			vector <int> depth(numberOfVertices);
			REP(u,1-1,numberOfVertices-1){
				depth[u] = 0;
			}
			
			DFS_Depth(depth, root, UNDEF);
			
			return depth;
		}
		
		
		vector <int> computeSubtreesSizes(int root){
			
			int UNDEF = -1;
			
			vector <int> subtreeSize(numberOfVertices);
			REP(u,1-1,numberOfVertices-1){
				subtreeSize[u] = 1;
			}
			
			DFS_SubtreeSize(subtreeSize, root, UNDEF);
			
			return subtreeSize;
		}
		
		
		LL ** FloydWarshall(){
			
			LL INF = (LL)1000000099 * 1000000009;
			
			LL **shortest = new LL *[numberOfVertices];
			REP(i,0,numberOfVertices-1){
				shortest[i] = new LL [numberOfVertices];
			}
			
			REP(u,0,numberOfVertices-1){
				REP(v,0,numberOfVertices-1){
					shortest[u][v] = INF;
				}
				
				shortest[u][u] = 0;
			}
			
			
			REP(i,0,numberOfEdges-1){
				int u = Edges[i].getDirectedEdge1().getStartVertex().getVertexID();    --u;
				int v = Edges[i].getDirectedEdge1().getEndVertex().getVertexID();    --v;
				int weight = Edges[i].getWeight();
				
				shortest[v][u] = shortest[u][v] = (   ( shortest[u][v] < weight ) ? shortest[u][v] : weight   );
			}
			
			REP(k,0,numberOfVertices-1){
				REP(u,0,numberOfVertices-1){
					REP(v,0,numberOfVertices-1){
						if( shortest[u][k] + shortest[k][v] < shortest[u][v] ){
							shortest[u][v] = shortest[u][k] + shortest[k][v];
						}
					}
				}
			}
			
			return shortest;
		}
		
		
		vector <LL> Dijkstra(int startVertex){
			
			const LL INF = (LL)1000000099 * 1000000009;
			
			--startVertex;
			
			PQPLL PQ;
			VLL shortest(numberOfVertices);
			bool visited[numberOfVertices];
			
			REP(u,0,numberOfVertices-1){
				shortest[u] = INF;
				visited[u] = false;
			}
			
			shortest[startVertex] = 0;
			PQ.push( make_pair( -shortest[startVertex], startVertex ) );
			
			int u, numberOfNeighbours, v, weight; 
			while( !PQ.empty() ){
				
				u = PQ.top().second;
				PQ.pop();
				
				if( visited[u] == false ){
					numberOfNeighbours = Graph[u].size();
				
					REP(i,0,numberOfNeighbours-1){
						v = Graph[u][i].getEndVertex().getVertexID();
						--v;
						
						weight = Graph[u][i].getWeight();
						
						if( shortest[u] + weight < shortest[v] ){
							shortest[v] = shortest[u] + weight;
							
							PQ.push( make_pair( -shortest[v], v ) );
						}	
					}
						
					visited[u] = true;
				}
			}
			
			return shortest;
		}
		
		
		/* 
			Min. 1 wierzcholek !!!
			Graf MUSI byc spojny !!! 
		*/
		vector <int> Kruskal(){
			
			#define MP make_pair
			
			
			vector <int> result;
			LL sum = 0;
			
			
			FindUnion findUnion(numberOfVertices);
			
			
			priority_queue <   pair<int, int>   > PQ;
			
			
			int weight, edgeID;
			
			REP(i,0,numberOfEdges-1){
				weight = Edges[i].getWeight();
				edgeID = Edges[i].getUndirectedEdgeID();
				
				PQ.push(   MP( -weight, edgeID )   );
			}
			
			
			int howManyEdgesChoosen = 0, u, v;
			
			while( howManyEdgesChoosen < numberOfVertices - 1 ){
				weight = -PQ.top().first;
				edgeID = PQ.top().second;
				
				u = Edges[edgeID-1].getDirectedEdge1().getStartVertex().getVertexID();    --u;
				v = Edges[edgeID-1].getDirectedEdge1().getEndVertex().getVertexID();    --v;
				
				PQ.pop();
				
				if( findUnion.Find(u) != findUnion.Find(v) ){
					findUnion.Union(u,v);
					
					result.push_back( edgeID );
					++howManyEdgesChoosen;
					
					sum += weight;
				}	
			}
			
			return result;
		}
		
		
		/* 
			Min. 1 wierzcholek !!!
			Graf MUSI byc spojny !!! 
		*/
		vector <int> Prim(int startVertex = 1){
			
			#define MP make_pair 
			
			const int INF = 1000000099;
			const int UNDEF = -1;
			
			vector <int> result;
			LL sum = 0;
			
			priority_queue <   pair< pair<int,int>, int >   > PQ; // (   (-waga_dojscia, wierzcholek), id_krawedzi   ) 
			
			bool isInTree[numberOfVertices];
			int leastWeight[numberOfVertices];
			
			REP(u,0,numberOfVertices-1){
				isInTree[u] = false;
				leastWeight[u] = INF;
			}
			leastWeight[startVertex-1] = 0;
			
			
			int u, numberOfNeighbours, v, weight, edgeID, currentNumberOfChoosenEdges = 0;
			
			
			PQ.push( MP( MP(-0, startVertex), UNDEF) );
			
			while( !PQ.empty()  &&  currentNumberOfChoosenEdges < numberOfVertices-1 ){
				weight = -PQ.top().first.first;
				u = PQ.top().first.second;    --u;
				edgeID = PQ.top().second;
				
				//cout << endl << endl << weight << " " << u << " " << edgeID << endl; 
				
				PQ.pop();
				
				if( !isInTree[u] ){
					//cout << "exploruje wierzcholek " << u << endl;
					
					isInTree[u] = true;
					
					if( u != startVertex-1 ){
						//cout << "dorzucam wierzcholek " << u << endl;
						result.push_back( edgeID );
						++currentNumberOfChoosenEdges;
						sum += weight;
					}
					/*else{
						cout << "nie dorzucam wierzcholka " << u << endl;
					}*/
					
					if( currentNumberOfChoosenEdges == numberOfVertices-1 ){
						//cout << "\nkoniec!\n";
						return result;
					}
					/*else{
						cout << "gramy dalej!\n";
					}*/
					
					
					numberOfNeighbours = Graph[u].size();
					//cout << numberOfNeighbours << endl;
					
					REP(i,0,numberOfNeighbours-1){
						v = Graph[u][i].getEndVertex().getVertexID();    --v;
						weight = Graph[u][i].getWeight();
						edgeID = Graph[u][i].getDirectedEdgeID();
						
						//cout << "sasiad nr " << i << ":\n" << weight << " " << v << " " << edgeID << endl; 
						
						if( !isInTree[v] ){
							//cout << "sasiad nr " << i << " nie jest jeszcze w drzewie!\n";
							
							//cout << leastWeight[v] << endl;
							
							if( weight < leastWeight[v] ){
								//cout << "wrzucam sasiada nr " << i << " do kolejki!\n";
								
								leastWeight[v] = weight;
								
								PQ.push(   MP( MP( -leastWeight[v], v+1 ), edgeID )   );
							}
							/*else{
								cout << "nie wrzucam sasiada nr " << i << " do kolejki!\n";
							}*/
						}
						/*else{
							cout << "sasiad nr " << i << " jest juz w drzewie!\n";
						}*/
						
						//cout << endl;
					}
				}
			}
		}
		
		
		void outputPrim(){
			vector <int> spanningTreePrim = (*this).Prim();
	
			sort( spanningTreePrim.begin(), spanningTreePrim.end() );
			
			cout << "\nAlgorytm Prima - id wybranych krawedzi:\n";
			
			
			REP(i,0,(this->numberOfVertices)-2){
				cout << spanningTreePrim[i] << " ";
			}
			cout << endl;
		}
		
		
		void output(){
			
			cout << endl << "Graph:\n" << numberOfVertices << " " << numberOfEdges << endl;
			
			/*REP(u,0,numberOfVertices-1){
				int numberOfNeighbours = Graph[u].size();
				
				REP(i,0,numberOfNeighbours-1){
					/*int edgeID = Graph[u][i].getEdgeID();
					int v = Graph[u][i].getEndVertex().getVertexID();
					int weight = Graph[u][i].getWeight();
					
					cout << edgeID << " " << u+1 << " " << v << " " << weight << endl;
					
					Graph[u][i].output();
				}
			}*/
			
			REP(i,0,numberOfEdges-1){
				Edges[i].output();
			}
		}
		
			
};


int main(){
	
	/* Creating an empty graph */
	WeightedUndirectedGraph G;
	
	
	/* Creating a graph */
	G.createWeightedUndirectedGraph();
	
	
	/* Coherentness */
	bool coherent = G.isCoherent();
	if( coherent == true ){
		cout << "\nSpojny!!!\n";
	}
	else{
		cout << "\nNiespojny!!!\n";
	}
	
	
	/* Dijkstra */
	int N = G.getNumberOfVertices();
	int startVertex = 1;
	vector <LL> shortest = G.Dijkstra(startVertex);
	cout << "\nShortest paths from vertex nr " << startVertex << ":" << endl;
	REP(v,0,N-1){
		cout << "shortest[" << v+1 << "] = ";
		
		//if( shortest[v] == INF ){
		//	cout << "INF" << endl;
		//}
		//else{
			cout << shortest[v] << endl;
		//}
	}
	
	
	/* Kruskal */
	vector <int> spanningTreeKruskal = G.Kruskal();
	
	sort( spanningTreeKruskal.begin(), spanningTreeKruskal.end() );
	
	cout << "\nAlgorytm Kruskala - id wybranych krawedzi:\n";
	
	int n = G.getNumberOfVertices();
	REP(i,0,n-2){
		cout << spanningTreeKruskal[i] << " ";
	}
	cout << endl;
	
	
	/* Prim */
	/*vector <int> spanningTreePrim = G.Prim();
	
	sort( spanningTreePrim.begin(), spanningTreePrim.end() );
	
	cout << "\nAlgorytm Prima - id wybranych krawedzi:\n";
	
	int nn = G.getNumberOfVertices();
	REP(i,0,nn-2){
		cout << spanningTreePrim[i] << " ";
	}
	cout << endl;*/
	
	G.outputPrim();
	
	
	/* Depths */
	int root = 1;
	vector <int> depth = G.computeDepths(root);
	int NN = G.getNumberOfVertices();
	cout << "\nDepths (root = " << root << "):\n";
	REP(u,0,NN-1){
		cout << "depth[" << u+1 << "] = " << depth[u] << endl;
	}
	cout << endl;
	
	
	/* Subtrees' sizes */
	int Root = 1;
	vector <int> subtreeSize = G.computeSubtreesSizes(Root);
	int NNN = G.getNumberOfVertices();
	cout << "\nSubtrees' sizes (root = " << Root << "):\n";
	REP(u,0,NNN-1){
		cout << "subtreeSize[" << u+1 << "] = " << subtreeSize[u] << endl;
	}
	cout << endl;
	
	
	/* Floyd-Warshall */
	LL **Shortest = G.FloydWarshall(), nnnn = G.getNumberOfVertices();
	
	cout << "\nFloyd-Warshall:\n";
	REP(u,0,nnnn-1){
		REP(v,0,nnnn-1){
			cout << Shortest[u][v] << " ";
		}
		cout << endl;
	}
	cout << endl;
	
	/* Output */
	G.output();
	
	
	// ------------------------------------
	
	
	
	/* Creating an empty graph */
	WeightedUndirectedGraph H;
	
	
	/* Creating a random graph */
	int maxNumberOfVertices = 6; 
	int maxNumberOfEdges = 10;
	int maxWeight = 10;
	
	srand( time(NULL) );
	H.createRandomWeightedUndirectedGraph(maxNumberOfVertices, maxNumberOfEdges, maxWeight);
	
	
	/* Floyd-Warshall */
	LL **best = H.FloydWarshall(), Hn = H.getNumberOfVertices();
	
	cout << "\nFloyd-Warshall:\n";
	REP(u,0,Hn-1){
		REP(v,0,Hn-1){
			cout << best[u][v] << " ";
		}
		cout << endl;
	}
	cout << endl;
	
	
	/* Output */
	H.output();
	
	return 0;
}
