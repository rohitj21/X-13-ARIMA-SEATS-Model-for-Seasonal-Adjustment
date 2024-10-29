#include<bits/stdc++.h>
#define mat(t) vector<vector<t>>
using namespace std;

string solve(int n,vector<int> &a){
    int i =1;
    int j=-1;
    int f=1;
    for(auto v:a){
        if(v==i){
            j++;
            i++;continue;
        }
        else{
            if(f){
                j=v;
                f=0;
                continue;
            }
            if(v!=j+1)return "NO";
            j++;
        }
    }

    return "YES";
}

int main(){
    vector<int> a={4,5,1,2,3};
    cout<<solve((int)a.size(),a)<<endl;
    return 0;
}

