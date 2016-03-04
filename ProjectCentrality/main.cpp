#include <iostream>
//#include"GraphMatrix.h"
#include"GraphMatrix.cpp"
using namespace std;

int main()
{
    Graph T;
    string name=".\\data\\1.txt";
    char *p=(char *)name.data();
    T.CreateGraphWithEdge(p);
   // T.Random_Walk_Closeness();
   T.Random_Walk_Betweeness();
   // T.LINE_RANK();
}
