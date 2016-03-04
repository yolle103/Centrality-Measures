#ifndef GRAPHMATRIX_H_INCLUDED
#define GRAPHMATRIX_H_INCLUDED
/*--------------------------------------------------------------------------------------
        Created by Liu jiamei on 2016/2/28
        Copyright 2016 Jilin University. All rights reserved.
        A measure for centrality
        including :


        1.betweenness centrality using random walk

            [1] M.E. J. Newman, A measure of betweenness centrality based on random walks,
            Social Networks, Volume 27, Issue 1, January 2005, Pages 39-54, ISSN 0378-8733,



        2.closeness centrality using random walk
            [2] Jae Dong N, Heiko R. Random walks on complex networks.
                [J]. Physical Review Letters, 2004, 92(11):87-93.



        3.betweenness centrality using line rank
            [3] Kang U, Papadimitriou S, Sun J, et al. Centralities in Large Networks: Algorithms and Observations[J].
                SDM, 2011:119-130.



        4.random walk centrality
            [4] https://github.com/alexrutherford/network_centrality





----------------------------------------------------------------------------------------*/


#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <ostream>
#include <stdio.h>
#include <sstream>
#include <time.h>
#include <string.h>
#include <cstring>
#include <ctime>
#include <vector>
#include <string>
#include <algorithm>
#include <list>
#include <map>
#include <ctype.h>
#include <ctime>
//#include"matrix.cpp"

#define inf 100000
using namespace std;


struct NodeinMatrix{
    float weight;
};

struct Adjacency_List
{
    string from,to;
};



class Graph
{

private:
    vector<vector<NodeinMatrix> > MyGraph;
    vector<vector<NodeinMatrix> > NGraph;
    map<string,float> Closeness_Centrality;
    map<string,float> Betweenness_Centrality;
    vector<Adjacency_List> List;
    map<string,int> Nodename;
    vector<float> Degree;
    vector<vector<NodeinMatrix> > SIM;
    vector<vector<NodeinMatrix> > TIM;
    vector<float> LineRank;
    vector<vector<NodeinMatrix> > Transposition(vector<vector<NodeinMatrix> > );
    bool Converge(vector<NodeinMatrix> &,vector<NodeinMatrix> &);
    void Load_Matrix();
    vector<NodeinMatrix>  Matrix_Mul(vector<vector<NodeinMatrix> > &, vector<NodeinMatrix>  &);
    vector<vector<NodeinMatrix> > Matrix_Mul_M(vector<vector<NodeinMatrix> > ,vector<vector<NodeinMatrix> > );
    vector<vector<NodeinMatrix> > Matrix_Add(vector<vector<NodeinMatrix> > &,vector<vector<NodeinMatrix> > &);
    int Graphsize;
    int Edgesize;
    void Inverse_Matrix(vector<vector<NodeinMatrix> > &);

    float Cal_Bet(int,vector<vector<NodeinMatrix> >);

public:
    void Random_walk();
    void CreateGraphWithEdge(char *);
    void Random_Walk_Closeness();
    void Random_Walk_Betweeness();
    void LINE_RANK();
    void PrintMatrix(vector<vector<NodeinMatrix> > a);
    Graph()
    {
        Graphsize=0;
    }
};


#endif // GRAPHMATRIX_H_INCLUDED
