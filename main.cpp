#include"libs/RootStructure.h"

/*
* This is a test file to ensure all headers are correctly set on your system.
*/
using namespace std;
using namespace Representation;
using namespace boost;
using namespace Eigen;

int main(int argc, char const ** argv)
{

	if(argc > 1){
		std::string loglevel = argv[1];
		Diagnositics::Log::ReportingLevel() = Diagnositics::Log::FromString(loglevel);
	}
	else{
		Diagnositics::Log::ReportingLevel() = Diagnositics::logINFO;

	}
    typedef Matrix<rational<int>, Dynamic, Dynamic> MatrixXr;



    LieBase<GroupType::E> f(6);

    // for (const auto& i : f.weightTower())
    // {
    //     std::cout<<i.omega.transpose()<<std::endl;
    // }

    typedef rational<int> frac;
    auto b = f.get_CoCartan();cout <<b << endl;
    auto a = PseudoInverse<frac>(f.get_Omega()).transpose();


    cout << a << endl;
    return 0;
}