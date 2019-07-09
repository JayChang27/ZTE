#define _CRT_SECURE_NO_WARNINGS

#include <iostream>
#include<cstdio>
#include<cstring>
#include<vector>
#include <stdlib.h>
#include<list>
#include<stdexcept>
#include<cassert>
using namespace std;
//存放原始拓扑信息
vector<int> gridTopo[956];
//存放所有业务需求信息
vector<int> request[4001];
//存放结果
vector<int> result[2001];
//业务数目
#define jobNum 1000
#define INF (~(0x1<<31))

// 构造“有向边”类
class DirectedEdge
{
public:
	DirectedEdge(int v = 0, int w = 0, int weight = 0, int width = 0)
	{
		this->v = v;
		this->w = w;
		this->weight = weight;
		this->width = width;
	}
	~DirectedEdge()
	{

	}
	int from()
	{
		return v;
	}
	int to()
	{
		return w;
	}
	int Weight()
	{
		return weight;
	}
	int Width()
	{
		return width;
	}
private:
	int v;
	int w;
	int weight;
	int width;
};

// 堆排序
class IndexMinPQ
{
public:
	IndexMinPQ(int maxN)
	{
		this->maxN = maxN;
		n = 0;
		keys = static_cast< vector<int> >(vector<int>(maxN + 1));
		pq = vector<int>(maxN + 1);
		qp = vector<int>(maxN + 1);
		for (int i = 0; i <= maxN; i++)
		{
			qp[i] = -1;
		}
	}
	~IndexMinPQ()
	{
		keys.~vector();
		pq.~vector();
		qp.~vector();
	}
	bool isEmpty()
	{
		return n == 0;
	}
	bool contains(int i)
	{
		return qp[i] != -1;
	}
	void insert(int i, int key)
	{
		n++;
		qp[i] = n;
		pq[n] = i;
		keys[i] = key;
		swim(n);
	}
	int delMin()
	{
		int min = pq[1];
		exch(1, n--);
		sink(1);
		qp[min] = -1;
		pq[n + 1] = -1;
		return min;
	}
	void decreaseKey(int i, int key)
	{
		keys[i] = key;
		swim(qp[i]);
	}
private:
	int maxN;
	int n;
	vector<int> pq;
	vector<int> qp;
	vector<int> keys;
	bool greater(int i, int j)
	{
		if (pq[i] > pq[j])
		{
			return true;
		}
		else
		{
			return false;
		}
	}
	void exch(int i, int j)
	{
		int tmp = pq[i];
		pq[i] = pq[j];
		pq[j] = tmp;
		qp[pq[i]] = i;
		qp[pq[j]] = j;
	}
	void swim(int k)
	{
		while (k > 1 && greater(k / 2, k))
		{
			exch(k, k / 2);
			k = k / 2;
		}
	}
	void sink(int k)
	{
		while (2 * k <= n)
		{
			int j = 2 * k;
			if (j < n && greater(j, j + 1))
			{
				j++;
			}
			if (!greater(k, j))
			{
				break;
			}
			exch(k, j);
			k = j;
		}
	}
};

// 构造“图”类
class Digraph
{
private:
	int V; // 顶点数
	int E; // 边数
	vector< vector<DirectedEdge> > adj; // 邻接矩阵
public:
	// 创建图
	Digraph(vector<int> gridTopo[])
	{
		V = gridTopo[0][0];
		int NE = gridTopo[0][1];
		E = 0;
		adj.resize(V);
		for (int i = 0; i < NE; i++)
		{
			int v = gridTopo[i + 1][0];
			int w = gridTopo[i + 1][1];
			int weight = gridTopo[i + 1][3];
			int width = gridTopo[i + 1][2] * 8 / 10;
			DirectedEdge e1(v, w, weight, width);
			DirectedEdge e2(w, v, weight, width);
			addEdge(e1);
			addEdge(e2);
		}
	}

	~Digraph()
	{
		adj.~vector();
	}
	void addEdge(DirectedEdge e)
	{
		int v = e.from();
		int w = e.to();
		adj[v].push_back(e);
		E++;
	}
	// 走过一个业务，这条边上带宽就减去走过业务的质量
	void decreaseEdge(DirectedEdge e, int w)
	{
		int s = e.from();
		int f = e.to();
		int width = e.Width() - w;
		int weight = e.Weight();
		DirectedEdge e_new1(s, f, weight, width);
		DirectedEdge e_new2(f, s, weight, width);
		int m = adj[s].size();
		for (int i = 0; i < m; i++)
		{
			int f_ = adj[s][i].to();
			if (f_ == f)
			{
				adj[s][i] = e_new1;
				break;
			}
		}
		int n = adj[f].size();
		for (int i = 0; i < n; i++)
		{
			int s_ = adj[f][i].to();
			if (s_ == s)
			{
				adj[f][i] = e_new2;
				break;
			}
		}
	}
	int getV() const
	{
		return V;
	}
	int getE() const
	{
		return E;
	}
	vector<DirectedEdge> Adj(int v)
	{
		return adj[v];
	}
};
// Dikstra算法类
class Dijkstra
{
public:
	Dijkstra(Digraph *G, int s, int width)
	{
		int n = G->getV();
		distTo_ = vector<int>(n);
		edgeTo = vector<DirectedEdge>(n);
		for (int v = 0; v < n; v++)
		{
			distTo_[v] = INF;
		}
		distTo_[s] = 0;
		pq = new IndexMinPQ(n);
		pq->insert(s, distTo_[s]);
		while (!pq->isEmpty())
		{
			int v = pq->delMin();
			for (DirectedEdge e : G->Adj(v))
			{
				relax(e, width);
			}
			
		}
		delete pq;
	}
	~Dijkstra()
	{
		delete pq;
	}
	int distTo(int v)
	{
		return distTo_[v];
	}
	bool hasPathTo(int v)
	{
		return distTo_[v] < INF;
	}
	list<DirectedEdge> pathTo(int v)
	{
		if (!hasPathTo(v))
		{
			return list<DirectedEdge>();
		}
		list<DirectedEdge> path;
		for (DirectedEdge e = edgeTo[v]; e.Weight() != 0; e = edgeTo[e.from()])
		{
			path.push_front(e);
		}
		return path;
	}
private:
	vector<int> distTo_; // distTo[v]表示起点s到终点v的最短距离
	vector<DirectedEdge> edgeTo; // edgeTo[v]表示s到v最短距离的最后一条边
	IndexMinPQ *pq; // 顶点的优先队列
	// 边的松弛
	void relax(DirectedEdge e, int width)
	{
		int v = e.from(), w = e.to();
		int wid = e.Width();
		if (distTo_[w] > distTo_[v] + e.Weight() && width <= wid)
		{
			distTo_[w] = distTo_[v] + e.Weight();
			edgeTo[w] = e;
			if (pq->contains(w))
			{
				pq->decreaseKey(w, distTo_[w]);
			}
			else
			{
				pq->insert(w, distTo_[w]);
			}
		}
	}
};

void clearVec()
{
	for (int i = 0; i < 956; i++) gridTopo[i].clear();
	for (int i = 0; i < 4001; i++) request[i].clear();
}
void readTxt()
{
	char readLine[1000];
	const char *delim = " ";
	char *p;
	for (int i = 0; i < 956; i++)
	{
		cin.getline(readLine, 1000);
		p = strtok(readLine, delim);
		while (p)
		{
			gridTopo[i].push_back(atoi(p));
			p = strtok(NULL, delim);
		}
	}
	for (int i = 0; i < 4001; i++)
	{
		cin.getline(readLine, 1000);
		p = strtok(readLine, delim);
		while (p)
		{
			request[i].push_back(atoi(p));
			p = strtok(NULL, delim);
		}
	}
}

int main()
{
	clearVec();
	//1.输入
	readTxt();

	result[0].push_back(0);
	Digraph *G = new Digraph(gridTopo);
	for (int i = 1; i <= jobNum; i++)
	{
		int width = request[4 * i - 3][1];
		int s = request[4 * i - 2].front();					//起点
		int f = request[4 * i - 2].back();					//终点
		Dijkstra *sp = new Dijkstra(G, s, width);
		if (sp->hasPathTo(f))
		{
			result[0][0] += sp->distTo(f)*width;
			result[2 * i - 1].push_back(i - 1);
			result[2 * i - 1].push_back(request[4 * i - 3][1]);
			result[2 * i].push_back(s);
			for (DirectedEdge e : sp->pathTo(f))
			{
				result[2 * i].push_back(e.to());
				G->decreaseEdge(e, width);
			}
		}
		else
		{
			result[0][0] = 0;
			break;
		}
		delete sp;
	}
	delete G;
	if (result[0][0] == 0)
	{
		cout << "NA";
	}
	else
	{
		for (int i = 0; i < 2001; i++)
		{
			int n = result[i].size();
			for (int j = 0; j < n; j++)
			{
				cout << result[i][j] << " ";
			}
			cout << endl;
		}
	}
	system("pause");
	return 0;
}
