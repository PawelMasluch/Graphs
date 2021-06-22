#include<iostream>
#include<vector>
#include<queue>
#include<new>


using namespace std;


typedef long long LL;

typedef vector <LL> VLL;

typedef pair <LL,LL> PLL;

typedef vector <PLL> VPLL;

typedef priority_queue < PLL > PQPLL;


#define REP(i,a,b) for(int i=a; i<=b; ++i)


const LL INF = (LL)1000000099 * 1000000009;


class Vertex{
	private:
		int vertexID;
		
	public:
		Vertex(int vertexID = 0){
			this -> vertexID = vertexID;
		}
		
		~Vertex() {}
		
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
};


class DirectedEdge{
	private:
		int directedEdgeID;
		Vertex startVertex;
		Vertex endVertex;
		int weight;
		
	public:
		DirectedEdge(int directedEdgeID, Vertex startVertex, Vertex endVertex, int weight = 1){
			this -> directedEdgeID = directedEdgeID;
			this -> startVertex = startVertex;
			this -> endVertex = endVertex;
			this -> weight = weight;
		}
		
		~DirectedEdge() {}
		
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


class WeightedDirectedGraph
{
	private:
		int numberOfVertices;
		int numberOfEdges;
		vector <DirectedEdge> *Graph;
		
		void operator += (DirectedEdge e){
			Graph[ e.getStartVertex().getVertexID() - 1 ].push_back( e );
		}
	
	public:
		
		WeightedDirectedGraph(){
			numberOfVertices = 0;
			numberOfEdges = 0;
			Graph = NULL;
		}
		
		
		~WeightedDirectedGraph(){
			
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
		
		
		void createWeightedDirectedGraph(){
			
			int numberOfVertices;
			cin >> numberOfVertices;
			
			int numberOfEdges;
			cin >> numberOfEdges;
			
			this -> numberOfVertices = numberOfVertices;
			this -> numberOfEdges = numberOfEdges;
			
			//cout << this -> numberOfVertices << " " << this -> numberOfEdges << endl;
			
			Graph = new vector <DirectedEdge> [this -> numberOfVertices];
			
			REP(directedEdgeID,1,this -> numberOfEdges){
				DirectedEdge e = DirectedEdge::createWeightedDirectedEdge(directedEdgeID);
				//cout << e.getStartVertex().getVertexID() << " " << e.getEndVertex().getVertexID() << " " << e.getWeight() << endl;
				*this += e;
			}
		}
		
		
		vector <LL> Dijkstra(int startVertex){
			
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
		
		
		void output(){
			
			cout << endl << "Graph:\n" << numberOfVertices << " " << numberOfEdges << endl;
			
			REP(u,0,numberOfVertices-1){
				int numberOfNeighbours = Graph[u].size();
				
				REP(i,0,numberOfNeighbours-1){
					/*int edgeID = Graph[u][i].getEdgeID();
					int v = Graph[u][i].getEndVertex().getVertexID();
					int weight = Graph[u][i].getWeight();
					
					cout << edgeID << " " << u+1 << " " << v << " " << weight << endl;
					*/
					Graph[u][i].output();
				}
			}	
		}
		
			
};


int main(){
	
	/* Creating an empty graph */
	WeightedDirectedGraph G;
	
	/* Creating a graph */
	G.createWeightedDirectedGraph();
	
	/* Dijkstra */
	int N = G.getNumberOfVertices();
	int startVertex = 1;
	vector <LL> shortest = G.Dijkstra(startVertex);
	cout << "\nShortest paths from vertex nr " << startVertex << ":" << endl;
	REP(v,0,N-1){
		cout << "shortest[" << v+1 << "] = ";
		
		if( shortest[v] == INF ){
			cout << "INF" << endl;
		}
		else{
			cout << shortest[v] << endl;	
		}
	}
	
	/* Output */
	G.output();
	
	return 0;
}
