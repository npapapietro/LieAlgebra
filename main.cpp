#include"libs/RootStructure.h"

/*
* This is a test file to ensure all headers are correctly set on your system.
*/
using namespace std;
using namespace Representation;
using namespace boost;

int main(int argc, char const ** argv)
{

	if(argc > 1){
		std::string loglevel = argv[1];
		Diagnositics::Log::ReportingLevel() = Diagnositics::Log::FromString(loglevel);
	}
	else{
		Diagnositics::Log::ReportingLevel() = Diagnositics::logINFO;

	}


    LieBase<GroupType::A> f(5);

    //cout << "Positive Roots\n";
    f.weightTower();

    // for (const auto& i: f.weightTower())
    // {
    //     cout<<int_cast_Vector(i.omega).transpose()<<endl;
    // }

    
    return 0;
}