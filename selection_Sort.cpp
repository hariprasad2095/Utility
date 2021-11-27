/******************************************************************************

                              Online C++ Compiler.
               Code, Compile, Run and Debug C++ program online.
Write your code in this editor and press "Run" button to compile and execute it.

*******************************************************************************/

#include <iostream>

using namespace std;
void selection_sort (int *,int);
int find_smallest (int *,int,int,int);
int
main ()
{
  int N;
  cout << "enter N\n";
  cin >>N;
  int *unsorted = new int (N);

  for (int i = 0; i < N; i++)
    {
      cout << "enter array value\n";
      cin >> *(unsorted + i);
    }
  selection_sort (unsorted,N);
  return 0;
}


void selection_sort(int* array,int N)
{
    int min_idx=0;
    for(int j=0;j<N;j++)
    {
        int temp=0;
        min_idx=j;
        min_idx = find_smallest(array,N,j,min_idx);
        cout<<"min_idx"<<min_idx;
        
        temp=*(array+min_idx);
        *(array+min_idx) = *(array+j);
        *(array+j) = temp;
    }
    cout <<"sorted array\n";
    for (int i = 0; i <N; i++)
    {
      cout <<*(array + i)<<"\n";
    }
    
}


int find_smallest(int* array, int N, int j,int min_idx)
{
    for(int i=j+1;i<N;i++)
    {
       if(*(array+i)<*(array+min_idx))
       {
           min_idx = i;
       }
    }
    return(min_idx);
} 
