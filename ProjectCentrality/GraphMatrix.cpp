
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
#include"GraphMatrix.h"

#define inf 100000
using namespace std;


void Graph::PrintMatrix(vector<vector<NodeinMatrix> > a)
{
    for(int i=0;i<a.size();i++)
    {
        for(int j=0;j<a[i].size();j++)
            cout<<a[i][j].weight<<"  ";
        cout<<endl;
    }
}

// multiply a matrix with a vector
vector<NodeinMatrix>  Graph::Matrix_Mul(vector<vector<NodeinMatrix> > &A,vector<NodeinMatrix>  &B)
{

    vector<NodeinMatrix>  R;
    R.resize(A.size());
    for(int i=0;i<A.size();i++)
    {
            R[i].weight=0;
            for(int x=0;x<B.size();x++)
            {
                R[i].weight+=A[i][x].weight*B[x].weight;
            }
    }
    return R;
}
vector<vector<NodeinMatrix> > Graph::Matrix_Mul_M(vector<vector<NodeinMatrix> > A,vector<vector<NodeinMatrix> > B)
{
    vector<vector<NodeinMatrix> > R;
    R.resize(A.size());
    for(int i=0;i<R.size();i++)
    {
        R[i].resize(B[0].size());
    }
    for(int i=0;i<A.size();i++)
    {
        for(int j=0;j<B[0].size();j++)
        {
             R[i][j].weight=0;
            for(int x=0;x<B.size();x++)
            {
                R[i][j].weight+=A[i][x].weight*B[x][j].weight;
            }
        }

    }
    return R;


}


//add two matrix
vector<vector<NodeinMatrix> > Graph::Matrix_Add(vector<vector<NodeinMatrix> > &A,vector<vector<NodeinMatrix> > &B)
{
    vector<vector<NodeinMatrix> > R;
    R=A;
    for(int i=0;i<A.size();i++)
    {
        for(int j=0;j<A[i].size();j++)
        {
            R[i][j].weight=A[i][j].weight+B[i][j].weight;
        }

    }
    return R;
}


//transposition a matrix
vector<vector<NodeinMatrix> > Graph::Transposition(vector<vector<NodeinMatrix> >  T)
{
    vector<vector<NodeinMatrix> > R;
    int row=T.size();
    int column=T[0].size();
    R.resize(column);
    for(int i=0;i<R.size();i++)
    {
        R[i].resize(row);
    }

    for(int i=0;i<row;i++)
    {
        for(int j=0;j<column;j++)
        {
            R[j][i].weight=T[i][j].weight;
        }
    }
    return R;
}


//test if vector converge
bool Graph::Converge(vector<NodeinMatrix>  &A,vector<NodeinMatrix> &B)
{
    float N=0.5;
    float sum=0;
    for(int i=0;i<A.size();i++)
    {
        sum+=fabs(A[i].weight-B[i].weight);
    }
    if(sum<N)
        return true;
    else
        return false;
}



//create graph with edge file
void Graph::CreateGraphWithEdge(char * filename)
{
    //open edge file
    fstream File;
    File.open(filename);

    //read edge info calculate graph size
    while(!File.eof())
    {
        Adjacency_List newEdge;
        File>>newEdge.from>>newEdge.to;
        List.push_back(newEdge);
      //  map<string,int>::iterator iter;
        if(Nodename.find(newEdge.from)==Nodename.end())
        {
            Nodename[newEdge.from]=Graphsize;
            Graphsize++;
        }
        if(Nodename.find(newEdge.to)==Nodename.end())
        {
            Nodename[newEdge.to]=Graphsize;
            Graphsize++;
        }
    }

    Edgesize=List.size();
    //convert to matrix form
   Load_Matrix();
}



void Graph::LINE_RANK()
{
    ofstream F;
    F.open(".\\result\\Linerank.txt");
    //initial sim and tim
    cout<<Edgesize<<endl;
   // float AA[Edgesize][Edgesize];
   // char a;
   // cin>>a;
    SIM.resize(Edgesize);
    TIM.resize(Edgesize);


    for(int i=0;i<Edgesize;i++)
    {
        SIM[i].resize(Graphsize);
        TIM[i].resize(Graphsize);
    }
 //   PrintMatrix(SIM);
   // PrintMatrix(TIM);

 cout<<"--------------------------------"<<endl;

    for(int i=0;i<Edgesize;i++)
    {
        SIM[i][Nodename[List[i].from]].weight=1;
        TIM[i][Nodename[List[i].to]].weight=1;
    }
  //

    //set the damping factor
    float c=0.85;

    //compute normalization factors
    vector<NodeinMatrix>  d1,d2,d,e;
    e.resize(SIM.size());

    for(int i=0;i<e.size();i++)
        e[i].weight=1;
    vector<vector<NodeinMatrix> > SIMT;
    SIMT=Transposition(SIM);
    d1=Matrix_Mul(SIMT,e);
    d2=Matrix_Mul(TIM,d1);
    d.resize(d2.size());
    for(int i=0;i<d2.size();i++)
    {
        if(d2[i].weight==0)
            d[i].weight=0;
        else
        {
            d[i].weight=(float)1/d2[i].weight;
        }

    }
    vector<NodeinMatrix> (d1).swap(d1);
    vector<NodeinMatrix> (d2).swap(d2);
    vector<NodeinMatrix> (e).swap(e);
    //random walk on TIMSIM'
  /*  cout<<"p1"<<endl;
    cout<<"d1"<<endl;
            for(int i=0;i<d1.size();i++)
            cout<<d1[i].weight<<endl;
            cout<<"d2"<<endl;
            for(int i=0;i<d2.size();i++)
            cout<<d2[i].weight<<endl;
            cout<<"d"<<endl;
            for(int i=0;i<d.size();i++)
            cout<<d[i].weight<<endl;*/


    vector<NodeinMatrix> v,v1,v2,v3,vc,LineRank;

    v.resize(Edgesize);

    vc.resize(Edgesize);
    //random initial vector v
    srand(time(NULL));
    for(int i=0;i<v.size();i++)
    {
        v[i].weight=rand()%10;
    }
    cout<<"p2"<<endl;
    vector<NodeinMatrix> (v1).swap(v1);
    vector<NodeinMatrix> (v2).swap(v2);
    vector<NodeinMatrix> (v3).swap(v3);

for(int i=0;i<v.size();i++)
            cout<<v[i].weight<<endl;


  //  while(!Converge(v,vc))
  for(int i=0;i<2;i++)
    {
      //  for(int i=0;i<v.size();i++)
        //    cout<<v[i].weight<<endl;
        cout<<"p3"<<endl;
        vc=v;
        v1=v;
        for(int i=0;i<v1.size();i++)
            v1[i].weight=d[i].weight*v[i].weight;

        v2=Matrix_Mul(SIMT,v1);
        v3=Matrix_Mul(TIM,v2);
        //add with the restart probability
        for(int i=0;i<v.size();i++)
        {
            v[i].weight=c*v3[i].weight+(float)1/Edgesize * (1-c);
        }
        vector<NodeinMatrix> (v1).swap(v1);
        vector<NodeinMatrix> (v2).swap(v2);
        vector<NodeinMatrix> (v3).swap(v3);

    }
    SIMT=Transposition(Matrix_Add(SIM,TIM));
    LineRank=Matrix_Mul(SIMT,v);
    for(int i=0;i<LineRank.size();i++)
        F<<LineRank[i].weight<<endl;
}



void Graph::Load_Matrix()
{
    //resize matrix with actual graph size
    MyGraph.resize(Graphsize);
    for( vector<vector<NodeinMatrix> >::iterator iter=MyGraph.begin();iter!=MyGraph.end();iter++)
    {
        (*iter).resize(Graphsize);
    }
    cout<<Graphsize<<endl;
    Degree.resize(Graphsize);
    //convert to matrix
    for(vector<Adjacency_List>::iterator iter=List.begin();iter!=List.end();iter++)
    {
        MyGraph[Nodename[(*iter).from]][Nodename[(*iter).to]].weight=1;
        Degree[Nodename[(*iter).from]]++;
        Degree[Nodename[(*iter).to]]++;
        MyGraph[Nodename[(*iter).to]][Nodename[(*iter).from]].weight=1;
    }

    //NGraph is a matrix containing the normalized factors
    NGraph.assign(MyGraph.begin(),MyGraph.end());
    for(int i=0;i<MyGraph.size();i++)
    {
        float sum=0;
        for(int j=0;j<NGraph.size();j++)
            sum+=NGraph[i][j].weight;
        for(int k=0;k<NGraph.size();k++)
            NGraph[i][k].weight/=sum;
    }
}





//calculate closeness centrality using random walk
void Graph::Random_Walk_Closeness()
{
    vector<NodeinMatrix> H;
    vector<vector<NodeinMatrix> > M;
    vector<vector<NodeinMatrix> > CM;

    CM.assign(MyGraph.begin(),MyGraph.end());

    for(int i=0;i<CM.size();i++)
    {
        for(int j=0;j<CM[i].size();j++)
        {
            if(i==j)
                CM[i][j].weight=1-CM[i][j].weight;
            else
                CM[i][j].weight=0-CM[i][j].weight;
        }
    }
cout<<"original graph"<<endl;
    for(int r=0;r<CM.size();r++)
    {
        for(int t=0;t<CM[r].size();t++)
            cout<<CM[r][t].weight<<"    ";
        cout<<endl;
        cout<<endl;
    }
    cout<<endl;

    ofstream F;
    F.open(".\\result\\closeness.txt");
    for(vector<vector<NodeinMatrix> >::iterator iter=CM.begin();iter!=CM.end();iter++)
    {
        //calculate the graph matrix delete its j-th row and column
        int row=iter-CM.begin();
        M.assign(CM.begin(),CM.end());
        M.erase(M.begin()+row);
        for(int cur=0;cur<M.size();cur++)
        {
            M[cur].erase(M[cur].begin()+row);
        }

        cout<<"delete "<<row<<"th row and column"<<endl;
         for(int r=0;r<M.size();r++)
        {
            for(int t=0;t<M[r].size();t++)
                cout<<M[r][t].weight<<" ";
            cout<<endl;
            cout<<endl;
        }

        Inverse_Matrix(M);
        cout<<endl;
        cout<<endl;
        cout<<"inverse"<<endl;
        for(int r=0;r<M.size();r++)
        {
            for(int t=0;t<M[r].size();t++)
                cout<<M[r][t].weight<<"  ";
            cout<<endl;
            cout<<endl;
        }
        cout<<endl;
        cout<<"-----------------------------------------------------------------------------------------------"<<endl;

        vector<float> H;
        for(int i=0;i<M.size();i++)
        {
            float sum=0;
            for(int j=0;j<M.size();j++)
                sum+=M[i][j].weight;
            H.push_back(sum);
        }
        float s=0;
        for(int k=0;k<H.size();k++)
        {
            s+=H[k];
        }

        F<<(float)(Graphsize)/s<<endl;
    }
}



//inverse matrix using LU method

void Graph::Inverse_Matrix(vector<vector<NodeinMatrix> > &A)
{
    vector<vector<NodeinMatrix> > L;
    vector<vector<NodeinMatrix> > U;
    vector<vector<NodeinMatrix> > r;
    vector<vector<NodeinMatrix> > u;

    int row,column;
    //initialing matrix L U r u
    L.resize(A.size());
    for( vector<vector<NodeinMatrix> >::iterator iter=L.begin();iter!=L.end();iter++)
    {
        (*iter).resize(A.size());
    }

    U.assign(L.begin(),L.end());

    r.assign(L.begin(),L.end());

    u.assign(L.begin(),L.end());


    for(int i=1;i<A.size();i++)
    {
         A[i][0].weight/=A[0][0].weight;
    }

  /*for(int i=0;i<A.size();i++)
    {
        for(int j=0;j<A.size();j++)
            cout<<A[i][j].weight<<"  ";
        cout<<endl;
    }*/


    float sum1,sum2;


    //calculate the L and U matrices
	for(int k=1;k<A.size();k++)
    {


        for(int j=k;j<A.size();j++)
        {
            sum1=0;
            for (int i=0;i<k;i++)
                sum1+=A[k][i].weight*A[i][j].weight;
            A[k][j].weight-=sum1;
		}


        for(int i=k+1;i<A.size();i++)
        {
            sum2=0;
            for(int j=0;j<k;j++)
                sum2+=A[i][j].weight*A[j][k].weight;
            A[i][k].weight=(A[i][k].weight-sum2)/A[k][k].weight;
        }


    }
        //Loading L and U matrices
		for(int i=0;i<A.size();i++)
			for(int j=0;j<A.size();j++)
			{
				if(i>j)
				{
					L[i][j].weight=A[i][j].weight;
					U[i][j].weight=0;
				}

				if(i==j)
				{
					U[i][j].weight=A[i][j].weight;
                    L[i][j].weight=1;
				}

				if(i<j)
                {
                    U[i][j].weight=A[i][j].weight;
                    L[i][j].weight=0;
                }

			}

				//inverse matrix U
				for (int i=0;i<A.size();i++)
				{
					u[i][i].weight=1/U[i][i].weight;
					for (int k=i-1;k>=0;k--)
					{
						sum1=0;
						for (int j=k+1;j<=i;j++)
							sum1+=U[k][j].weight*u[j][i].weight;
						u[k][i].weight=-sum1/U[k][k].weight;
					}
				}

				//inverse matrix L
				for (int i=0;i<A.size();i++)
				{
					r[i][i].weight=1;
					for (int k=i+1;k<A.size();k++)
					{
						for (int j=i;j<=k-1;j++)
							r[k][i].weight-=L[k][j].weight*r[j][i].weight;
					}
				}

				//r multiplied by u to produce the inverse matrix of A

				for(int i=0;i<A.size();i++)
				{
					for(int j=0;j<A.size();j++)
                        A[i][j].weight=0;
				}

				for(int i=0;i<A.size();i++)
				{
					for(int j=0;j<A.size();j++)
					{
						for(int k=0;k<A.size();k++)
						{
							A[i][j].weight+=u[i][k].weight*r[k][j].weight;
						}
					}
				}

 }

void Graph::Random_Walk_Betweeness()
{

    ofstream F;
    F.open(".\\result\\betweenness.txt");
   vector<vector<NodeinMatrix> > BGraph;
   BGraph.assign(MyGraph.begin(),MyGraph.end());
 //  cout<<"Degree"<<endl;
//   cout<<Degree[0]<<endl;
 //  for(int k=0;k<Degree.size();k++)
   // cout<<Degree[k]<<endl;
    //cout<<"why?-----------------------------"<<endl;
   //D-A
   for(int i=0;i<BGraph.size();i++)
   {
        for(int j=0;j<BGraph.size();j++)
        {
            if(i==j)
            {
                BGraph[i][j].weight=Degree[i]-BGraph[i][j].weight;
            }
            else
            {
                BGraph[i][j].weight=-BGraph[i][j].weight;
            }

        }
   }
   vector<vector<NodeinMatrix> > T;
   //remove row and column
    for(vector<vector<NodeinMatrix> >::iterator iter=BGraph.begin();iter!=BGraph.end();iter++)
    {

        int row=iter-BGraph.begin();
        T.assign(BGraph.begin(),BGraph.end());
        T.erase(T.begin()+row);
        for(int cur=0;cur<T.size();cur++)
        {
            T[cur].erase(T[cur].begin()+row);
        }

   /*    cout<<"delete "<<row<<"th row and column"<<endl;
         for(int r=0;r<T.size();r++)
        {
            for(int t=0;t<T[r].size();t++)
                cout<<T[r][t].weight<<" ";
            cout<<endl;
            cout<<endl;
        }
*/
        Inverse_Matrix(T);
  /*      cout<<endl;
        cout<<endl;
        cout<<"inverse"<<endl;
        for(int r=0;r<T.size();r++)
        {
            for(int t=0;t<T[r].size();t++)
                cout<<T[r][t].weight<<"  ";
            cout<<endl;
            cout<<endl;
        }
        cout<<endl;
        cout<<"-----------------------------------------------------------------------------------------------"<<endl;*/
       vector<NodeinMatrix>  Zeros;
        Zeros.resize(Graphsize);
        //back to normal
        NodeinMatrix M;
        M.weight=0;
        for(int cur=0;cur<T.size();cur++)
        {
            T[cur].insert(T[cur].begin()+row,M);
        }
     //  cout<<"back to normal"<<endl;
        T.insert(row+T.begin(),Zeros);
      /*   for(int r=0;r<T.size();r++)
        {
            for(int t=0;t<T[r].size();t++)
                cout<<T[r][t].weight<<"  ";
            cout<<endl;
            cout<<endl;
        }
        cout<<endl;*/
        //calculate betweenness
        vector<float> Betweenness;
        Betweenness.resize(Graphsize);
        Betweenness[row] = Cal_Bet(row,T);
        cout<<"betweenness:"<<endl;
        F<<Betweenness[row]<<endl;
      //  char k;
        //cin>>k;

    }
}


float Graph::Cal_Bet(int i,vector<vector<NodeinMatrix> > T)
{
    float result=0;
    float I_sum=0;
    float I_sum_sum=0;


    for(int s=0;s<Graphsize;s++)
    {

        for(int t=s+1;t<Graphsize;t++)
        {
            float I=0;
            //calculate I_sum
          //  cout<<"s  ="<<s<<"  t = "<<t<<endl;
            for(int j=0;j<Graphsize;j++)
            {
                if(i==s||i==t)
                {
                    I=1;
                    break;
                }
                else
                {
                    I+=(MyGraph[i][j].weight)*fabs(T[i][s].weight-T[i][t].weight-T[j][s].weight+T[j][t].weight);
             //       cout<<"when j= "<<j<<" I ="<<I<<endl;
              //      cout<<fabs(T[i][s].weight-T[i][t].weight-T[j][s].weight+T[j][t].weight)<<endl;
               //     cout<<MyGraph[i][j].weight<<endl;
                }
            }

       //     cout<<"when s="<<s<<"  t= "<<t<<" we have I="<<I<<endl;
            I_sum+=I;

        }

        I_sum_sum+=I_sum;


    }
    result=I_sum_sum/(Graphsize*(Graphsize-1));
  //  cout<<"result for i "<<result<<endl;
    return result;
}


void Graph::Random_walk()
{
    //save it
    ofstream F;
    F.open(".\\result\\randomwalk.txt");

    int start=0;
    vector<float> H;
    H.resize(Graphsize);
    srand(time(NULL));
    //cast a random walk on the graph matrix and calculating time at every vertex
    for(int i=0;i<inf;i++)
    {
        H[start]++;
        int position=rand()%1000;
        float sum=0;
        int newPosition;
        float newsum=0;
        for(int i=0;i<Graphsize;i++)
        {
            newsum=sum;
            newsum+=MyGraph[start][i].weight*1000;
            if(sum<position&&newsum>=position&&MyGraph[start][i].weight!=0)
            {
                newPosition=i;
                break;
            }
            else
            {
                sum=newsum;

            }


        }
        start=newPosition;
    }
    for(int i=0;i<H.size();i++)
    {
        H[i]/=inf;
        H[i]*=100;
        F<<H[i]<<endl;
    }
}
