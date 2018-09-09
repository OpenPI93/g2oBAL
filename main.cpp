#include <iostream>
#include "balproblem.h"

using namespace std;

int main()
{
    balProblem myBAL("./data1.txt");
    myBAL.buildProblem();
    myBAL.solveProblem(40);
//    myBAL.writeToPLY("./data.ply");
    myBAL.showByPangolin();

    return 0;
}
