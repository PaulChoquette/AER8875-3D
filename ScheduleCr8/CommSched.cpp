
#include <fstream>
#include <iostream>
#include <string>
using namespace std;

string FileName = "Communication_Schedule.txt";
// int nzone = 7;
// int ConctMax = 4;
// int OffLoad[7][4] = {{2,5,-1,-1},
//                     {5,3,6,-1},
//                     {0,3,4,-1},
//                     {2,1,6,-1},
//                     {2,6,-1,-1},
//                     {0,1,6,-1},
//                     {1,3,4,5}};

int nzone = 4;
int ConctMax = 4;
int OffLoad[4][4] = {{1,3,2,-1},
                    {0,2,-1,-1},
                    {3,1,0,-1},
                    {2,0,-1,-1}};


int main() 
{
    int ** Zone2Zone;
    // Create input matrix
    
    Zone2Zone = new int*[nzone];
    for (int i = 0;i<nzone;++i) {
        Zone2Zone[i] = new int[nzone];
        for (int j = 0;j<ConctMax;++j) {
            Zone2Zone[i][j] = OffLoad[i][j];
        }
    }    

    // Real Code
    // Allocate Variables
    bool AnyTx = 1;
    int Rnd = 0;
    bool * isBusy;
    isBusy = new bool[nzone];
    int Zone2;
    ofstream file;
    file.open(FileName);

    while (AnyTx) {
        
        // Reset is buzy
        for (int i=0;i<nzone;++i) {
            isBusy[i] = 0;
        }
        AnyTx = 0;

        for (int i=0;i<nzone;++i) {
            for (int j=0;j<ConctMax;++j) {
                Zone2 = Zone2Zone[i][j];
                if ((0==isBusy[i])&&(Zone2!=-1)&&(0==isBusy[Zone2])) {
                    if (!AnyTx) {
                        file<<"RND "<<Rnd<<endl;
                        AnyTx = 1;
                    }
                    isBusy[i] = 1;
                    isBusy[Zone2] = 1;
                    Zone2Zone[i][j] = -1;
                    for (int j2=0;j2<ConctMax;++j2){
                        if (Zone2Zone[Zone2][j2]==i) {
                            Zone2Zone[Zone2][j2]=-1;
                            break;
                        }
                    }
                    // Output to file
                    file<<i<<" "<<Zone2<<endl;
                    break;
                }
            }
        }
        Rnd += 1;
    }
    file.close();




    // Delete allocated matrix 
    for (int i = 0;i<nzone;++i) {
        delete[] Zone2Zone[i];
    }
    delete[] Zone2Zone;
    delete[] isBusy;

return 1;}